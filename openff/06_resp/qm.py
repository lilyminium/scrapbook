
from qcfractal import interface as ptl

from info import PROTOCOLS, PCM_KEYWORDS


def run_geometry_optimization(client, qcmols, method="PW6B95", basis="cc-pV(D+d)Z", **kwargs):
    spec = {
        "keywords": None,
        "qc_spec": {
            "driver": "gradient",
            "method": method,
            "basis": basis,
            "program": "psi4",
            "protocols": PROTOCOLS
        },
    }

    # Ask the server to compute a new computation
    response = client.add_procedure("optimization", "geometric", spec, qcmols)
    return response.ids


def compute_energy(client, qcmols, method="PW6B95", basis="cc-pV(D+d)Z", state="gas"):
    keywords = {"maxiter": 300}
    if state != "gas":
            keywords.update(PCM_KEYWORDS)

    keyword_id = client.add_keywords([ptl.models.KeywordSet(values=keywords)])[0]

    computation = dict(
        program="psi4",
        basis=basis,
        method=method,
        driver="energy",
        molecule=qcmols,
        keywords=keyword_id,
        protocols=PROTOCOLS,
    )
    response = client.add_compute(**computation)
    return response.ids