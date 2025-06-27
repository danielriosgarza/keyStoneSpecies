# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 10:35:41 2025

@author: drgarza
"""

import os
# Create directory for the license
os.makedirs("licenses", exist_ok=True)

def create_gurobi_license():
    license_content = (
        "# Gurobi WLS license file\n"
        "# Your credentials are private and should not be shared or copied to public repositories.\n"
        "# Visit https://license.gurobi.com/manager/doc/overview for more information.\n"
        "WLSACCESSID=d2dc56cd-66c9-4422-a1e1-9fe2d819fbb5\n"
        "WLSSECRET=c6b325ef-a1d8-4c4c-a4f4-6d5e25ac1f34\n"
        "LICENSEID=940603"
    )
    with open("licenses/gurobi.lic", "w") as f:
        f.write(license_content)
    print("License file created at licenses/gurobi.lic")

if __name__ == "__main__":
    create_gurobi_license()