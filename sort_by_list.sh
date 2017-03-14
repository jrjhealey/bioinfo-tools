#!/bin/bash

#!/bin/bash
FILE_TO_SORT="$1"
INDEX_FILE="$2"
TMP_FILE=$(mktemp)

while read LINE; do
    grep "$LINE" "$FILE_TO_SORT" >>"$TMP_FILE"
done <"$INDEX_FILE"

mv -f "$TMP_FILE" "$FILE_TO_SORT"
