import re
def sync_reference_name(ref_name):
    ref_name = ref_name.upper()
    match = re.search("(?<=CHR).*", ref_name)
    if match:
        ref_name = match.group(0)
    if ref_name == "M":
        ref_name = "MT"
    if '_' in ref_name:
        ref_name = ''
    return ref_name