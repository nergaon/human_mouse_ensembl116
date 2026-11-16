"""Supplementary table numbers of the six-dataset version: the coverage table (internal S7) was removed, so internal
Tables S8-S13 are published as S7-S12. Internal names are kept in the scripts; renum() converts text and file names
when they are written."""
import re

OLD2NEW = {n: n for n in range(1, 7)}
OLD2NEW.update({n: n - 1 for n in range(8, 14)})
REMOVED = {7}


def new_number(n):
    n = int(n)
    if n in REMOVED:
        return None
    return OLD2NEW.get(n, n)


def name(table_name):
    """'Table_S8' -> 'Table_S7'; 'S8' -> 'S7'; None for the removed table."""
    m = re.fullmatch(r'(Table_)?S(\d+)', table_name)
    if not m:
        return table_name
    k = new_number(m.group(2))
    return None if k is None else f'{m.group(1) or ""}S{k}'


LIST = re.compile(r'(Tables? S\d+(?:(?:, | and |-|–)S\d+)*|Table_S\d+)')


def renum(text):
    if not isinstance(text, str):
        return text

    def fix(m):
        def one(k):
            v = new_number(k.group(1))
            if v is None:
                raise ValueError('reference to the removed Table S7: ' + text[max(0, m.start() - 80):m.end() + 20])
            return f'S{v}'
        return re.sub(r'S(\d+)', one, m.group(0))
    return LIST.sub(fix, text)
