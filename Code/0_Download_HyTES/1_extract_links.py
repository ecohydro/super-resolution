### This file extracts the email contents from the apple script and writes the zip file links to a text file (text.txt)
#!/usr/bin/env python3
import sys
import re

# Sets contents of email to email_body variable
email_body = sys.argv[1]

# Directory where we wish to store the text file
dir = '/Users/rstofer/documents/Github/super-resolution/Data/'

# only searches for .zip links
link = 'https?://.*.zip'

# Searches for urls
urls = re.findall(link, email_body)

# Write to text file
with open(dir + 'text.txt','a') as f:
    [f.write(url) for url in urls]
    f.write("\n")
