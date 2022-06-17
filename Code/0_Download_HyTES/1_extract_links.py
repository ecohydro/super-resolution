#!/usr/bin/env python3
import sys
import re

email_body = sys.argv[1]

dir = '/Users/rstofer/documents/Github/super-resolution/Data/'

link = 'https?://.*.zip'

urls = re.findall(link, email_body)

with open(dir + 'text.txt','a') as f:
    [f.write(url) for url in urls]
    f.write("\n")
