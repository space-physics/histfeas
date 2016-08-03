try:
    from pathlib import Path
    Path().expanduser()
else:
    from pathlib2 import Path
