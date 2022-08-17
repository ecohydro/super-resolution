#!/usr/bin/env python3
import sys
import re

email_body = sys.argv[1]

dir = '/Users/rrstofer/documents/Github/super-resolution/Outputs'

urls = re.findall('http?://(?:[-\w.]|(?:%[\da-fA-F]{2}))+', email_body)

with open(dir + 'text.txt','w') as f:
    [f.write(url) for url in urls]