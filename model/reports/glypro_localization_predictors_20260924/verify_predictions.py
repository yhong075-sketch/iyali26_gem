"""Validate official downloads; never substitute a pending job with a prediction."""
from pathlib import Path
import csv,json,hashlib,zipfile,datetime
ROOT=Path(__file__).resolve().parent
q="".join(x.strip() for x in (ROOT/"target.fasta").read_text().splitlines() if not x.startswith(">"))
SHA="d7a158663b2c1d752f1b14bec8cdf30f135470992eb8478b26521964cb159b59"
assert len(q)==454 and hashlib.sha256(q.encode()).hexdigest()==SHA
identity="YALI1E16433g_W29_AOW05368.1"
rows=list(csv.DictReader((ROOT/"deeploc/results_6AB5CD1E00184CB90257F738.csv").open()))
assert len(rows)==1 and rows[0]["Protein_ID"]==identity
attention=list(csv.DictReader((ROOT/"deeploc/alpha_yali1e16433g_w29_aow053681.csv").open()))
assert "".join(x["AA"] for x in attention)==q
dl_json=json.loads((ROOT/"deeploc/results.json").read_text())
assert dl_json["info"]=={"size":1,"failedjobs":0} and list(dl_json["sequences"])==[identity]
thresholds=dict(zip(dl_json["Localization"]+dl_json["Membrane_types"],map(float,dl_json["Threshold"]+dl_json["Threshold_memtype"])))
assert len(thresholds)==14
scores={k:float(rows[0][k]) for k in thresholds}
assert {k for k in list(thresholds)[:10] if scores[k]>thresholds[k]}==set(rows[0]["Localizations"].split("|"))
assert {k for k in list(thresholds)[10:] if scores[k]>thresholds[k]}==set(rows[0]["Membrane types"].split("|"))
result={"query_sha256":SHA,"verified_at":datetime.datetime.now(datetime.timezone.utc).isoformat(),"deeploc":{"version":"2.1","job_id":"6AB5CD1E00184CB90257F738","status":"completed_verified","mode":"High-quality (Slow), ProtT5; long output","localizations":rows[0]["Localizations"].split("|"),"membrane_types":rows[0]["Membrane types"].split("|"),"signals":rows[0]["Signals"].split("|"),"scores_csv":scores,"thresholds":thresholds,"threshold_source":"Official results.json loaded by the job page; archived after browser resource inventory discovery; thresholds parsed directly","returned_sequence_prefix_length_verified":454}}
for tool,name,n in [("signalp","prediction_results.txt",70),("targetp","output_protein_type.txt",200)]:
 j=json.loads((ROOT/tool/"output.json").read_text())
 assert j["INFO"]["failedjobs"]==0 and j["INFO"]["size"]==1 and list(j["SEQUENCES"])==[identity]
 raw=[line.split("\t") for line in (ROOT/tool/name).read_text().splitlines() if line and not line.startswith("#")]
 assert len(raw)==1 and raw[0][0]==identity and raw[0][1]=="OTHER" and raw[0][-1]==""
 with zipfile.ZipFile(ROOT/tool/"output_all_results.zip") as z:
  assert z.testzip() is None
  names=[v for v in z.namelist() if v.endswith(".txt") and "_YALI" in v]
  assert len(names)==1
  text=z.read(names[0]).decode(); pos=[r.split("\t") for r in text.splitlines() if r and not r.startswith("#")]
  assert len(pos)==n and "".join(r[1] for r in pos)==q[:n] and [int(r[0]) for r in pos]==list(range(1,n+1))
  (ROOT/tool/names[0]).write_text(text)
 result[tool]={"version":"6.0" if tool=="signalp" else "2.0","job_id":"6AB5CD3500184CEE6E75D672" if tool=="signalp" else "6AB5CD4100184D1721C360D1","status":"completed_verified","organism":j["ORG"],"mode":j.get("MODE"),"prediction":raw[0][1],"scores":dict(zip(["OTHER","SP(Sec/SPI)"] if tool=="signalp" else ["OTHER","SP","mTP"],map(float,raw[0][2:-1]))),"cleavage_site":None,"returned_sequence_prefix_length_verified":n,"full_input_length":454,"limitation":"Input was full pinned sequence; downloadable per-residue output covers the stated N-terminal prefix, not full sequence."}
assert result["signalp"]["mode"]=="slow-sequential" and result["signalp"]["organism"]=="eukarya"
assert result["targetp"]["organism"]=="Non-Plant"
tm=ROOT/"deeptmhmm"
meta=json.loads((tm/".biolib/metadata.json").read_text())
assert meta["exit_code"]==0
inputs=list((tm/"biolib-input-files").glob("*.fasta")); assert len(inputs)==1
rawinput="".join(x.strip() for x in inputs[0].read_text().splitlines() if not x.startswith(">"))
assert rawinput==q
lines=[x.strip() for x in (tm/"predicted_topologies.3line").read_text().splitlines() if x.strip()]
assert len(lines)==3 and lines[0]==">"+identity+" | GLOB" and lines[1]==q and lines[2]=="I"*454
assert "Number of predicted TMRs: 0" in (tm/"TMRs.gff3").read_text()
probs=list(csv.reader((tm/(identity+"_probs.csv")).read_text().splitlines()[2:]))
assert len(probs)==454 and "".join(r[0].split()[1] for r in probs)==q
result["deeptmhmm"]={"version":"1.0.57","job_id":"1736a09d-6406-41fc-8c79-695754acd5a5","status":"completed_verified","prediction":"GLOB","predicted_transmembrane_regions":0,"topology":"I at residues 1-454","returned_sequence_prefix_length_verified":454,"input_file_full_sequence_verified":True,"official_metadata":meta,"limitation":"GLOB/inside is a topology classification; not direct evidence of native cytoplasm vs nucleus, organelle residency, or catalytic orientation."}
(ROOT/"prediction_verification.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps(result,indent=2))
