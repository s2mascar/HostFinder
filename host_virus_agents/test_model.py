from transformers import AutoTokenizer, AutoModelForCausalLM
import torch
import time

MODEL_PATH = "./models/Qwen3-0.6B"

print("Loading tokenizer...")
tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH)

print("Loading model...")
model = AutoModelForCausalLM.from_pretrained(
    MODEL_PATH,
    dtype=torch.float32
)

print("Model loaded.")

messages = [
    {
        "role": "system",
        "content": """
You are a scientific host-virus analysis agent.

Your job is to examine a host and virus pair.

Return ONLY valid JSON in this format:

{
    "host": "...",
    "virus": "...",
    "possible_interaction": true,
    "reason": "..."
}
"""
    },
    {
        "role": "user",
        "content": """
Host: Tribolium castaneum
Virus: Hubei partiti-like virus 31
"""
    }
]

text = tokenizer.apply_chat_template(
    messages,
    tokenize=False,
    add_generation_prompt=True,
    enable_thinking=False
)

inputs = tokenizer(text, return_tensors="pt")

print("Generating...")

start = time.time()

with torch.no_grad():
    outputs = model.generate(
        **inputs,
        max_new_tokens=150,
        do_sample=False
    )

print(f"Generation took {time.time() - start:.1f} seconds")

response = tokenizer.decode(
    outputs[0][inputs.input_ids.shape[1]:],
    skip_special_tokens=True
)

print("\nMODEL RESPONSE:")
print(response)