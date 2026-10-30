import sys, json, urllib.request, urllib.parse
q = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 6
url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q, "format": "json", "resultType": "core", "pageSize": n})
d = json.load(urllib.request.urlopen(url, timeout=60))
for r in d["resultList"]["result"]:
    print("PMID", r.get("pmid"), "PMC", r.get("pmcid"), r.get("pubYear"), r.get("journalInfo", {}).get("journal", {}).get("title"), "DOI", r.get("doi"))
    print("  T:", r.get("title")); print("  A:", r.get("authorString", "")[:120])
    print("  AB:", (r.get("abstractText") or "")[:1800].replace("<h4>", " ").replace("</h4>", " ")); print()
