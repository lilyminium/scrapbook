import os


array_values = []
for i in range(1, 4067):
    if not os.path.isfile(f"remappings/mol{i:04d}.smi"):
    #     continue
    # if not os.path.isfile(f"charges/at_am1bcc/mol{i:04d}.csv"):
        array_values.append(str(i))

print(",".join((array_values)))