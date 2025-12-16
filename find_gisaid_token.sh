#!/bin/bash

TOKEN_PATH=$(python -c "import outbreak_data, os; print(os.path.join(os.path.dirname(os.path.dirname(outbreak_data.__file__)), '.Python_outbreak_info_token.txt'))" 2>/dev/null)

if [ -z "$TOKEN_PATH" ]; then
    echo "Error: Could not determine token file path. Is outbreak_data installed?" >&2
    exit 1
fi

if [ ! -f "$TOKEN_PATH" ]; then
    echo "Error: Token file does not exist at: $TOKEN_PATH" >&2
    echo "Please run 'python gisaid_authentication.py' to authenticate first." >&2
    exit 1
fi

echo "$TOKEN_PATH"

