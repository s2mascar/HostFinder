<#
Export an immutable historical benchmark snapshot without importing the pipeline.
Run from the project root. No model, network, or classification changes.
Null means not recorded; historical cache records are context, not proof of use.
#>
[CmdletBinding()]
param(
    [string]$InputCsv = 'test_pairs.csv',
    [string]$ResultsCsv = 'results/bioresearch_env_v06_results.csv',
    [string]$ContextCache = 'biological_context_cache_v06.json',
    [string]$VirusCache = 'virus_taxonomy_cache_v06.json',
    [string]$OutputJson = 'results/benchmark_diagnostic_v06.json'
)
$ErrorActionPreference = 'Stop'

function Read-JsonField($Row, [string]$Name) {
    $property = $Row.PSObject.Properties[$Name]
    if ($null -eq $property -or [string]::IsNullOrWhiteSpace($property.Value)) { return $null }
    return ($property.Value | ConvertFrom-Json)
}
function Normalize-Name([string]$Name) {
    return (($Name.ToLowerInvariant() -replace '[^a-z0-9]+', ' ').Trim() -replace '\s+', ' ')
}
function Cache-Entry($Cache, [string]$Key) {
    if ($null -eq $Cache) { return $null }
    $property = $Cache.PSObject.Properties[$Key]
    if ($null -ne $property) { return $property.Value }
    return $null
}
function Read-OptionalCache([string]$Path) {
    if (Test-Path -LiteralPath $Path -PathType Leaf) {
        return (Get-Content -LiteralPath $Path -Raw | ConvertFrom-Json)
    }
    return $null
}

# Do not overwrite historical results, inputs, or previous diagnostic runs.
if (Test-Path -LiteralPath $OutputJson) { throw "Output already exists: $OutputJson" }
$definitions = @(Import-Csv -LiteralPath $InputCsv)
$results = @(Import-Csv -LiteralPath $ResultsCsv)
if ($definitions.Count -ne $results.Count) { throw 'Input/result row counts differ.' }
$contexts = Read-OptionalCache $ContextCache
$viruses = Read-OptionalCache $VirusCache
$records = @()
for ($index = 0; $index -lt $definitions.Count; $index++) {
    $definition = $definitions[$index]
    $row = $results[$index]
    foreach ($field in @('host', 'virus', 'expected_status')) {
        if ([string]::IsNullOrWhiteSpace($definition.$field) -or $definition.$field -cne $row.$field) {
            throw "Definition/result mismatch at row $($index + 1), field $field"
        }
    }
    $correct = $row.predicted_status -ceq $definition.expected_status
    if ([bool]::Parse($row.correct) -ne $correct) { throw "Incorrect stored correctness flag at row $($index + 1)" }
    $papers = @(Read-JsonField $row 'paper_diagnostics')
    if ($papers.Count -eq 1 -and $null -eq $papers[0]) { $papers = @() }
    if ($papers.Count -ne [int]$row.papers_analyzed) { throw "Missing paper diagnostics at row $($index + 1)" }
    if (@($papers | Group-Object paper_id | Where-Object Count -gt 1).Count) { throw 'Duplicate episode paper ID.' }
    $exact = @($papers | Where-Object classification -eq 'EXACT_SUPPORT' | ForEach-Object paper_id)
    $related = @($papers | Where-Object classification -eq 'TARGET_HOST_RELATED' | ForEach-Object paper_id)
    if ($exact.Count -ne [int]$row.exact_supporting_papers -or $related.Count -ne [int]$row.related_papers) {
        throw "Evidence count mismatch at row $($index + 1)"
    }
    $contextKey = (Normalize-Name $row.host) + '||' + (Normalize-Name $row.virus)
    $context = Cache-Entry $contexts $contextKey
    $virusKey = 'name:' + (Normalize-Name $row.virus)
    if ($null -ne $context -and $context.target_virus.tax_id) { $virusKey = 'taxid:' + $context.target_virus.tax_id }
    $paperRecords = @()
    foreach ($paper in $papers) {
        $paperRecords += [ordered]@{
            paper_id = $paper.paper_id
            pmid = $paper.pmid
            pmcid = $paper.pmcid
            title = $paper.title
            classification = $paper.classification
            classification_basis = $paper.classification_basis
            extraction_status = $paper.extraction_status
            extraction_error = $paper.extraction_error
            evidence_type = [ordered]@{
                host_relationship_type = $paper.host_edge.relationship_type
                comparison_relationship_type = $paper.comparison_edge.relationship_type
                provenance = 'model labels and verifier diagnostics, not independent annotation'
            }
            evidence_scope = [ordered]@{
                support_mode = $paper.host_edge.support_mode
                natural_vs_experimental = $null
                status = 'not separately recorded; inspect relationship types and supporting text'
            }
            target_host_binding = $paper.host_edge
            target_virus_binding = [ordered]@{
                reported_name = $paper.host_virus_name
                matches_target = $paper.host_edge.target_virus_match
                aliases_used = $paper.target_virus_aliases
                tax_id = $paper.target_virus_tax_id
                scientific_name = $paper.target_virus_scientific_name
            }
            comparison_binding = $paper.comparison_edge
            original_diagnostic = $paper
        }
    }
    $records += [ordered]@{
        input_row = $index + 1
        host = $row.host
        virus = $row.virus
        expected_status = $definition.expected_status
        expected_evidence = $definition.expected_evidence
        test_type = $definition.test_type
        reference_hint = $definition.reference_hint
        predicted_status = $row.predicted_status
        correct = $correct
        result_origin = 'historical_artifact; no live run or extraction replay'
        taxonomy_resolution = [ordered]@{
            historical_context_snapshot = $context
            historical_virus_snapshot = (Cache-Entry $viruses $virusKey)
            status = 'snapshot context only; original run cache identity is not recorded'
        }
        host_aliases_used = $null
        virus_aliases_used = @($papers | ForEach-Object {
            [ordered]@{ paper_id = $_.paper_id; aliases = $_.target_virus_aliases }
        })
        literature_queries_attempted = $null
        search_status = [ordered]@{
            query_count_recorded = [int]$row.searches
            per_source_attempts = $null
            complete = $null
            reason = 'query strings, HTTP outcomes and source-success counts absent from CSV'
        }
        papers_retrieved = [int]$row.papers_retrieved
        papers_retrieved_semantics = 'selected candidates, not raw retrieval count or total hits'
        raw_retrieved_count = $null
        candidates_before_selection = $null
        candidate_papers_considered = $paperRecords
        exact_support_papers = $exact
        related_only_papers = $related
        confidence = $row.confidence
        final_reason = $row.reason
        runtime_seconds = [double]::Parse($row.runtime_seconds, [cultureinfo]::InvariantCulture)
        cache_usage = [ordered]@{ hits = $null; misses = $null; status = 'not recorded' }
        original_result = $row
    }
}
$paths = @($InputCsv, $ResultsCsv, $ContextCache, $VirusCache, 'AGENTS.md', 'PROJECT_SPEC.md',
    'ARCHITECTURE_REVIEW.md', 'evaluate_pairs.py', 'evaluate_env.py', 'search_agent.py',
    'literature_search.py', 'taxonomy_aliases.py', 'evidence_agent.py', 'judge_agent.py',
    'requirements.txt', 'run_benchmark.sh', 'run_benchmark_v2.sh', 'export_benchmark_diagnostics.ps1')
$paths += @(Get-ChildItem bioresearch_env -Filter '*.py' -File | ForEach-Object { 'bioresearch_env/' + $_.Name })
$manifest = @($paths | Sort-Object -Unique | ForEach-Object {
    $hash = $null
    if (Test-Path -LiteralPath $_ -PathType Leaf) { $hash = (Get-FileHash -LiteralPath $_ -Algorithm SHA256).Hash.ToLowerInvariant() }
    [ordered]@{ path = $_; sha256 = $hash }
})
$payload = [ordered]@{
    schema_version = 1
    mode = 'historical_artifact_export'
    reproducibility = 'deterministic export; current source hashes do not authenticate historical executable/model versions'
    unavailable = @('raw papers', 'model inputs/outputs', 'host aliases used', 'query strings',
        'per-source success/failure', 'raw candidate pool', 'original model revision', 'cache-hit events')
    manifest = $manifest
    summary = [ordered]@{
        total = $records.Count
        correct = @($records | Where-Object { $_.correct }).Count
        failed_extractions = @($records | ForEach-Object { $_.candidate_papers_considered } | Where-Object extraction_status -eq 'FAILED').Count
    }
    cases = $records
}
# CreateNew refuses overwrites even if a competing invocation creates the path.
$stream = [IO.File]::Open($OutputJson, [IO.FileMode]::CreateNew, [IO.FileAccess]::Write)
try {
    $text = ($payload | ConvertTo-Json -Depth 60) + "`n"
    $bytes = [Text.UTF8Encoding]::new($false).GetBytes($text)
    $stream.Write($bytes, 0, $bytes.Length)
} finally { $stream.Dispose() }
Write-Output "Exported $($records.Count) historical cases; $($payload.summary.correct)/$($records.Count) correct -> $OutputJson"
