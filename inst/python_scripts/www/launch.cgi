#!/usr/bin/python
import matplotlib
# before anything else can load, force matplotlib to not use any
# Xwindows backend.
matplotlib.use('Agg')

import sys
import cgi
import os
import tempfile
import SuperMASSA

print 'Content-type: text/html\r\n\r\n'

# redirect all errors to stdout
sys.stderr = sys.stdout

try:
    sys.path.append('/usr/lib/cgi-bin/')
    import inference
except Exception as e:
    print e

form = cgi.FieldStorage()

options_list = []
for f in form:
    out = form[f].value
    if f != 'file' and f != 'f1_parent_data':
        if out != 'off':
            options_list.append( '--' + str(f) )
            options_list.append( form[f].value )
        else:
            options_list.append( '--' + str(f) )

# write the file to a temporary location
uploadedfile = form['file']
temp, temp_filename = tempfile.mkstemp(text = True)
temp_output = os.fdopen(temp, 'w')
file_data = uploadedfile.file.read()
if len( file_data ) != 0:
    temp_output.write(file_data)
    options_list.extend(['--file', temp_filename])
temp_output.close()

# if the parent data is included, write that to a temporary location
if form['f1_parent_data'] != '':
    uploadedfile = form['f1_parent_data']
    temp_parental, temp_parental_filename = tempfile.mkstemp(text = True)
    temp_output = os.fdopen(temp_parental, 'w')
    parent_file_data = uploadedfile.file.read()
    if parent_file_data != '':
        # empty file = no file
        temp_output.write(parent_file_data)
        temp_output.close()
        options_list.extend(['--f1_parent_data', temp_parental_filename])

SuperMASSA.cgi_main(options_list)
