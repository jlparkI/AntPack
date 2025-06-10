"""Sets up license keys for antpack."""
import os


def run_license_key_setter():
    """Asks the user for a license key and email and
    writes it to an appropriate file."""
    print("Enter your license key (see the docs for instructions "
            "to obtain one):")
    license_key = input()
    print("Enter your email:")
    email = input()

    fpath = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(fpath, "license_keys.txt"), "w+",
            encoding="utf-8") as fhandle:
        fhandle.write("\t".join([license_key, email]))

    print("License key has been saved.")


def get_license_key_info():
    """Retrieves a saved license key."""
    fpath = os.path.abspath(os.path.dirname(__file__))
    if not os.path.exists(os.path.join(fpath, "license_keys.txt")):
        raise RuntimeError("License key not found! Please run AntPack-setup "
                "first. If you need a license key, please see the docs.")
    with open(os.path.join(fpath, "license_keys.txt"), "r",
            encoding="utf-8") as fhandle:
        license_key, email = fhandle.readline().strip().split()

    return license_key, email
