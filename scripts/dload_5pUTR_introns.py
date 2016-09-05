#!/usr/bin/env python
from intermine.webservice import Service
service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")
query = service.new_query("FivePrimeUTRIntron")
query.add_view("orf.secondaryIdentifier")

for row in query.rows():
    print row["orf.secondaryIdentifier"]
