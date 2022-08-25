### This file extracts the email contents from the apple script and writes the zip file links to a text file (links.txt)
#!/usr/bin/env python3
import sys
import re

# Sets contents of email to email_body variable
email_body = sys.argv[1]

# Directory where we wish to store the text file
dir = '/Users/rstofer/documents/Github/super-resolution/Data/Raw'

# only searches for .zip links
link = 'https?://.*.zip'

# Searches for urls
urls = re.findall(link, email_body)

# Write to text file
with open(dir + 'links.txt','a') as f:
    [f.write(url) for url in urls]
    f.write("\n")
