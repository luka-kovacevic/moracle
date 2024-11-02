from transformers import AutoTokenizer, AutoModelForCausalLM
import peft
import torch
import os 

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)

access_token = os.environ.get("HF_TOKEN") # have your huggingface token saved as an environment variable $HF_TOKEN
model_id = "meta-llama/Llama-2-7b-chat-hf"
lora_id = "moleculargpt_lora"

model = AutoModelForCausalLM.from_pretrained(model_id).to(device)
tokenizer = AutoTokenizer.from_pretrained(model_id)

model.load_adapter(lora_id)

print(model.device)



# breakpoint()
#### need to change the stuff below to conform to the moleculargpt prompt formatting

messages = [
    {"role": "system", "content": "You are a pirate chatbot who always responds in pirate speak!"},
    {"role": "user", "content": "Who are you?"},
]

input_ids = tokenizer.apply_chat_template(
    messages,
    add_generation_prompt=True,
    return_tensors="pt"
).to(model.device)

terminators = [
    tokenizer.eos_token_id,
    tokenizer.convert_tokens_to_ids("<|eot_id|>")
]

outputs = model.generate(
    input_ids,
    max_new_tokens=256,
    eos_token_id=terminators,
    do_sample=True,
    temperature=0.6,
    top_p=0.9,
)
response = outputs[0][input_ids.shape[-1]:]
print(tokenizer.decode(response, skip_special_tokens=True))
