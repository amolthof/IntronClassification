import sys
import re

def parse_line(line):
    # Use regex to find key-value pairs
    pattern = re.compile(r'(\w+) "([^"]*)"')
    matches = pattern.findall(line)

    attributes = {k: v for k, v in matches}

    return attributes

def main():
    # Read lines from stdin
    lines = [line.strip() for line in sys.stdin if line.strip()]

    # Parse each line
    entries = [parse_line(line) for line in lines]

    # Get all unique keys in the order they first appear
    all_keys = []
    for entry in entries:
        for key in entry:
            if key not in all_keys:
                all_keys.append(key)

    # Print headers with tabs
    print('\t'.join(all_keys))

    # Print each entry with tabs
    for entry in entries:
        values = [entry.get(key, '') for key in all_keys]
        print('\t'.join(values))

if __name__ == "__main__":
    main()
