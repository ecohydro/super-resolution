#!/usr/bin/env python3

## This file extracts the url links from each e-mail and writes it to a text file

import sys
import re

email_body = sys.argv[1]

dir = '/Users/rstofer/documents/Github/super-resolution/Data/'

link = 'https?://.*.zip'

urls = re.findall(link, email_body)

with open(dir + 'text.txt','a') as f:
    [f.write(url) for url in urls]
    f.write("\n")
