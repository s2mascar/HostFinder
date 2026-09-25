import os

os.environ.setdefault(
    "PYTORCH_CUDA_ALLOC_CONF",
    "expandable_segments:True"
)

import gc
import json
import re
import sys
import time

import requests
import torch

from transformers import (
    AutoTokenizer,
    AutoModelForCausalLM
)

from literature_search import (
    search_pubmed,
    search_pmc,
    search_europe_pmc,
    merge_papers
)

from taxonomy_aliases import (
    get_host_aliases
)

from bioresearch_env.relationship_language import (
    interaction_proximity_score
)


# ============================================================
# SETTINGS
# ============================================================

MODEL_PATH = "./models/Qwen3-8B"

MAX_SEARCH_QUERIES = 3
MAX_RESULTS_PER_SOURCE = 5
MAX_CANDIDATE_PAPERS = 8

MAX_RETRIES = 4

NCBI_PAUSE_SECONDS = 1.0


# ============================================================
# MODEL
# ============================================================

DEVICE = (
    "cuda"
    if torch.cuda.is_available()
    else "cpu"
)

if DEVICE == "cuda":

    if torch.cuda.is_bf16_supported():
        DTYPE = torch.bfloat16
    else:
        DTYPE = torch.float16

    if hasattr(
        torch.backends.cuda,
        "enable_cudnn_sdp"
    ):
        torch.backends.cuda.enable_cudnn_sdp(
            False
        )

else:

    DTYPE = torch.float32


print(f"Loading model on: {DEVICE}")
print(f"Using dtype: {DTYPE}")

tokenizer = AutoTokenizer.from_pretrained(
    MODEL_PATH
)

model = AutoModelForCausalLM.from_pretrained(
    MODEL_PATH,
    dtype=DTYPE,
    low_cpu_mem_usage=True
)

model = model.to(
    DEVICE
)

model.eval()

print("Model loaded.")


# ============================================================
# GPU CLEANUP
# ============================================================

def cleanup_gpu():

    gc.collect()

    if torch.cuda.is_available():
        torch.cuda.empty_cache()


# ============================================================
# GENERATION
# ============================================================

def generate_text(
    messages,
    max_new_tokens=250,
    max_input_tokens=2048
):

    text = tokenizer.apply_chat_template(
        messages,
        tokenize=False,
        add_generation_prompt=True,
        enable_thinking=False
    )

    inputs = tokenizer(
        text,
        return_tensors="pt",
        truncation=True,
        max_length=max_input_tokens
    ).to(
        model.device
    )

    outputs = None

    try:

        with torch.inference_mode():

            outputs = model.generate(
                **inputs,
                max_new_tokens=max_new_tokens,
                do_sample=False,
                pad_token_id=tokenizer.eos_token_id
            )

        input_length = (
            inputs[
                "input_ids"
            ].shape[1]
        )

        generated = (
            outputs[0][input_length:]
            .detach()
            .cpu()
        )

        response = tokenizer.decode(
            generated,
            skip_special_tokens=True
        )

        del generated

        return response

    finally:

        if outputs is not None:
            del outputs

        del inputs

        cleanup_gpu()


# ============================================================
# JSON EXTRACTION
# ============================================================

def extract_json_object(
    text
):

    decoder = json.JSONDecoder()

    for i, char in enumerate(
        text
    ):

        if char != "{":
            continue

        try:

            obj, _ = decoder.raw_decode(
                text[i:]
            )

            if isinstance(
                obj,
                dict
            ):
                return obj

        except json.JSONDecodeError:
            continue

    raise ValueError(
        "No valid JSON object found."
    )


# ============================================================
# NAME NORMALIZATION
# ============================================================

def normalize_name(
    text
):

    if not text:
        return ""

    text = str(
        text
    ).lower()

    text = re.sub(
        r"[^a-z0-9]+",
        " ",
        text
    )

    return " ".join(
        text.split()
    )


def contains_name(
    text,
    name
):

    if not text or not name:
        return False

    return (
        normalize_name(
            name
        )
        in normalize_name(
            text
        )
    )


def contains_any_alias(
    text,
    aliases
):

    return any(
        contains_name(
            text,
            alias
        )
        for alias
        in aliases
    )


# ============================================================
# SEARCH ALIAS SELECTION
# ============================================================

def alias_priority(
    alias,
    host
):

    alias_norm = normalize_name(
        alias
    )

    host_norm = normalize_name(
        host
    )

    if alias_norm == host_norm:
        return -100

    # Skip abbreviations such as L. noctiluca.
    if re.match(
        r"^[A-Za-z]\.\s",
        alias
    ):
        return -50

    score = 0

    # Common names from NCBI are often lowercase.
    if alias and alias[0].islower():
        score += 20

    # Multi-word common names tend to be useful.
    if len(
        alias.split()
    ) >= 2:
        score += 5

    # Avoid extremely short aliases.
    if len(alias) < 5:
        score -= 20

    return score


def get_useful_aliases(
    host,
    host_aliases
):

    candidates = []

    for alias in host_aliases:

        if (
            normalize_name(alias)
            == normalize_name(host)
        ):
            continue

        if re.match(
            r"^[A-Za-z]\.\s",
            alias
        ):
            continue

        candidates.append(
            alias
        )

    candidates = sorted(
        candidates,
        key=lambda alias:
            alias_priority(
                alias,
                host
            ),
        reverse=True
    )

    return candidates


# ============================================================
# VIRUS BROADENING
# ============================================================

def shorten_virus_name(
    virus
):

    parts = virus.split()

    if (
        len(parts) > 1
        and parts[-1].isdigit()
    ):
        return " ".join(
            parts[:-1]
        )

    return virus


# ============================================================
# MODEL SEARCH PLANNING
# ============================================================

def generate_model_queries(
    host,
    virus,
    host_aliases
):

    aliases_text = "\n".join(
        f"- {alias}"
        for alias
        in host_aliases
    )

    system_prompt = """
You are a scientific literature search-planning agent.

Your only task is to suggest search terms.

You are NOT deciding whether an interaction exists.

Rules:

1. Never invent host synonyms.
2. Only use the validated host names supplied below.
3. Preserve the target virus name unless cautious broadening
   is useful.
4. Never infer that the host and virus interact.
5. Return no more than 3 searches.

Return ONLY valid JSON:

{
    "queries": [
        {
            "host_term": "...",
            "virus_term": "...",
            "reason": "..."
        }
    ]
}
"""

    user_prompt = f"""
TARGET HOST:
{host}

VALIDATED HOST NAMES:
{aliases_text}

TARGET VIRUS:
{virus}

Generate literature searches.
"""

    messages = [
        {
            "role": "system",
            "content":
                system_prompt
        },
        {
            "role": "user",
            "content":
                user_prompt
        }
    ]

    try:

        response = generate_text(
            messages,
            max_new_tokens=250,
            max_input_tokens=1600
        )

        data = extract_json_object(
            response
        )

        return data.get(
            "queries",
            []
        )

    except Exception as error:

        print(
            "Search query generation "
            f"failed: {error}"
        )

        return []


# ============================================================
# BUILD FINAL SEARCH PLAN
# ============================================================

def generate_search_queries(
    host,
    virus,
    host_aliases
):

    queries = []

    seen = set()

    def add_query(
        host_term,
        virus_term,
        reason
    ):

        if (
            len(queries)
            >= MAX_SEARCH_QUERIES
        ):
            return

        if not host_term:
            return

        if not virus_term:
            return

        key = (
            normalize_name(
                host_term
            ),
            normalize_name(
                virus_term
            )
        )

        if key in seen:
            return

        seen.add(
            key
        )

        queries.append({
            "host_term":
                host_term,

            "virus_term":
                virus_term,

            "reason":
                reason
        })

    # --------------------------------------------------------
    # Query 1: exact scientific name
    # --------------------------------------------------------

    add_query(
        host,
        virus,
        "Exact host and exact virus."
    )

    # --------------------------------------------------------
    # Query 2: best NCBI host alias
    # --------------------------------------------------------

    useful_aliases = (
        get_useful_aliases(
            host,
            host_aliases
        )
    )

    if useful_aliases:

        add_query(
            useful_aliases[0],
            virus,
            (
                "Validated NCBI host alias "
                "with the exact virus."
            )
        )

    # --------------------------------------------------------
    # Let Qwen fill remaining space.
    # --------------------------------------------------------

    model_queries = (
        generate_model_queries(
            host,
            virus,
            host_aliases
        )
    )

    valid_alias_norms = {
        normalize_name(alias)
        for alias
        in host_aliases
    }

    valid_alias_norms.add(
        normalize_name(
            host
        )
    )

    for item in model_queries:

        if (
            len(queries)
            >= MAX_SEARCH_QUERIES
        ):
            break

        host_term = str(
            item.get(
                "host_term",
                ""
            )
        ).strip()

        virus_term = str(
            item.get(
                "virus_term",
                ""
            )
        ).strip()

        if (
            normalize_name(
                host_term
            )
            not in valid_alias_norms
        ):
            continue

        add_query(
            host_term,
            virus_term,
            item.get(
                "reason",
                "Model-generated search."
            )
        )

    # --------------------------------------------------------
    # Deterministic fallback
    # --------------------------------------------------------

    if (
        len(queries)
        < MAX_SEARCH_QUERIES
    ):

        shortened = (
            shorten_virus_name(
                virus
            )
        )

        if shortened != virus:

            add_query(
                host,
                shortened,
                (
                    "Exact host with cautious "
                    "virus-name broadening."
                )
            )

    if (
        len(queries)
        < MAX_SEARCH_QUERIES
        and len(useful_aliases) > 1
    ):

        add_query(
            useful_aliases[1],
            virus,
            (
                "Second validated NCBI "
                "host alias."
            )
        )

    return queries


# ============================================================
# SOURCE RETRIES
# ============================================================

def call_source_with_retry(
    source_name,
    function,
    host_term,
    virus_term
):

    last_error = None

    for attempt in range(
        MAX_RETRIES
    ):

        try:

            results = function(
                host_term,
                virus_term,
                max_results=
                    MAX_RESULTS_PER_SOURCE
            )

            return (
                results,
                None
            )

        except requests.HTTPError as error:

            last_error = error

            status = (
                error.response.status_code
                if error.response
                is not None
                else None
            )

            if status not in {
                429,
                500,
                502,
                503,
                504
            }:
                break

            wait = 2 ** attempt

            print(
                f"{source_name} HTTP "
                f"{status}. Retrying "
                f"in {wait}s..."
            )

            time.sleep(
                wait
            )

        except requests.RequestException as error:

            last_error = error

            wait = 2 ** attempt

            time.sleep(
                wait
            )

        except Exception as error:

            last_error = error
            break

    return (
        [],
        str(last_error)
    )


# ============================================================
# SEARCH DATABASES
# ============================================================

def search_query_sources(
    host_term,
    virus_term
):

    all_results = []

    statuses = []

    sources = [
        (
            "PubMed",
            search_pubmed,
            True
        ),
        (
            "PMC",
            search_pmc,
            True
        ),
        (
            "Europe PMC",
            search_europe_pmc,
            False
        )
    ]

    for (
        source_name,
        function,
        is_ncbi
    ) in sources:

        results, error = (
            call_source_with_retry(
                source_name,
                function,
                host_term,
                virus_term
            )
        )

        if error is None:

            statuses.append({
                "source":
                    source_name,

                "success":
                    True,

                "error":
                    None
            })

            for paper in results:

                paper.setdefault(
                    "search_queries",
                    []
                )

                paper[
                    "search_queries"
                ].append({
                    "host_term":
                        host_term,

                    "virus_term":
                        virus_term
                })

            all_results.extend(
                results
            )

        else:

            statuses.append({
                "source":
                    source_name,

                "success":
                    False,

                "error":
                    error
            })

            print(
                f"WARNING: {source_name} "
                f"failed: {error}"
            )

        if is_ncbi:

            time.sleep(
                NCBI_PAUSE_SECONDS
            )

    return (
        all_results,
        statuses
    )


# ============================================================
# PAPER RANKING
# ============================================================

def candidate_score(
    paper,
    host_aliases,
    virus,
    host=None
):

    title = (
        paper.get(
            "title",
            ""
        )
        or ""
    )

    abstract = (
        paper.get(
            "abstract",
            ""
        )
        or ""
    )

    full_text = (
        paper.get(
            "full_text",
            ""
        )
        or ""
    )

    score = 0

    # ========================================================
    # TARGET VIRUS
    # ========================================================

    virus_title = contains_name(
        title,
        virus
    )

    virus_abstract = contains_name(
        abstract,
        virus
    )

    virus_full = contains_name(
        full_text,
        virus
    )

    # ========================================================
    # EXACT SCIENTIFIC HOST NAME
    # ========================================================

    exact_host_title = (
        contains_name(
            title,
            host
        )
        if host
        else False
    )

    exact_host_abstract = (
        contains_name(
            abstract,
            host
        )
        if host
        else False
    )

    exact_host_full = (
        contains_name(
            full_text,
            host
        )
        if host
        else False
    )

    # ========================================================
    # VALIDATED HOST ALIASES
    # ========================================================

    alias_title = contains_any_alias(
        title,
        host_aliases
    )

    alias_abstract = contains_any_alias(
        abstract,
        host_aliases
    )

    alias_full = contains_any_alias(
        full_text,
        host_aliases
    )

    # Virus evidence.
    if virus_title:
        score += 15

    if virus_abstract:
        score += 9

    if virus_full:
        score += 4

    # Exact scientific host names receive much more weight than
    # generic aliases such as "mouse" or "human".
    if exact_host_title:
        score += 15
    elif alias_title:
        score += 5

    if exact_host_abstract:
        score += 10
    elif alias_abstract:
        score += 4

    if exact_host_full:
        score += 4
    elif alias_full:
        score += 2

    # Strong co-occurrence bonuses.
    if (
        virus_title
        and exact_host_title
    ):
        score += 25

    if (
        virus_abstract
        and exact_host_abstract
    ):
        score += 18

    if (
        virus_abstract
        and alias_abstract
        and not exact_host_abstract
    ):
        score += 6

    # Most important new signal: host + virus + relationship
    # language in the same local passage.
    score += interaction_proximity_score(
        abstract,
        host,
        host_aliases,
        virus
    )

    # Full text is much longer and noisier, so use half-weight.
    score += (
        interaction_proximity_score(
            full_text,
            host,
            host_aliases,
            virus
        )
        // 2
    )

    if full_text:
        score += 2

    return score


# ============================================================
# RUN SEARCH AGENT
# ============================================================

def run_search_agent(
    host,
    virus
):

    host_aliases = (
        get_host_aliases(
            host
        )
    )

    queries = (
        generate_search_queries(
            host,
            virus,
            host_aliases
        )
    )

    print(
        "\n"
        + "=" * 80
    )

    print(
        "SEARCH PLAN"
    )

    print(
        "=" * 80
    )

    for i, query in enumerate(
        queries,
        start=1
    ):

        print(
            f'{i}. '
            f'"{query["host_term"]}" '
            f'AND '
            f'"{query["virus_term"]}"'
        )

        print(
            "   Reason:",
            query[
                "reason"
            ]
        )

    all_papers = []

    metadata = {
        "queries_planned":
            len(queries),

        "queries_executed":
            0,

        "source_successes":
            0,

        "source_failures":
            [],

        "host_aliases":
            host_aliases
    }

    # --------------------------------------------------------
    # IMPORTANT:
    # Execute ALL planned searches before ranking.
    #
    # Do not stop just because the first noisy query
    # returned >= 8 papers.
    # --------------------------------------------------------

    for i, query in enumerate(
        queries,
        start=1
    ):

        print(
            "\n"
            + "#" * 80
        )

        print(
            f"EXECUTING SEARCH "
            f"{i}/{len(queries)}"
        )

        print(
            "#" * 80
        )

        papers, statuses = (
            search_query_sources(
                query[
                    "host_term"
                ],
                query[
                    "virus_term"
                ]
            )
        )

        metadata[
            "queries_executed"
        ] += 1

        for status in statuses:

            if status[
                "success"
            ]:

                metadata[
                    "source_successes"
                ] += 1

            else:

                metadata[
                    "source_failures"
                ].append({
                    "query":
                        query,

                    "source":
                        status[
                            "source"
                        ],

                    "error":
                        status[
                            "error"
                        ]
                })

        all_papers.extend(
            papers
        )

    unique_papers = merge_papers(
        all_papers
    )

    unique_papers = sorted(
        unique_papers,
        key=lambda paper:
            candidate_score(
                paper,
                host_aliases,
                virus,
                host=host
            ),
        reverse=True
    )

    unique_papers = (
        unique_papers[
            :MAX_CANDIDATE_PAPERS
        ]
    )

    metadata[
        "candidate_papers"
    ] = len(
        unique_papers
    )

    metadata[
        "retrieval_complete"
    ] = (
        len(
            metadata[
                "source_failures"
            ]
        )
        == 0
    )

    print(
        "\n"
        + "=" * 80
    )

    print(
        "FINAL SEARCH-AGENT RESULTS"
    )

    print(
        "=" * 80
    )

    print(
        "Search strategies executed:",
        metadata[
            "queries_executed"
        ]
    )

    print(
        "Unique candidate papers retained:",
        len(
            unique_papers
        )
    )

    print(
        "Literature source failures:",
        len(
            metadata[
                "source_failures"
            ]
        )
    )

    return {
        "papers":
            unique_papers,

        "search_metadata":
            metadata,

        "host_aliases":
            host_aliases
    }


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    if len(sys.argv) != 3:

        print(
            'Usage: python search_agent.py '
            '"HOST" "VIRUS"'
        )

        sys.exit(1)

    result = run_search_agent(
        sys.argv[1],
        sys.argv[2]
    )

    for i, paper in enumerate(
        result[
            "papers"
        ],
        start=1
    ):

        print(
            "\n"
            + "=" * 80
        )

        print(
            f"CANDIDATE PAPER {i}"
        )

        print(
            "=" * 80
        )

        print(
            "Title:",
            paper.get(
                "title"
            )
        )

        print(
            "PMID:",
            paper.get(
                "pmid"
            )
        )

        print(
            "PMCID:",
            paper.get(
                "pmcid"
            )
        )