# extracts date from dardarfile name. Also available as a class method, 
# but havent figured out to use it without creating a class instance

import numpy as np
from datetime import datetime
import os

def filename2date(filename):
        filename = os.path.basename(filename)
        filename = filename.split("_")[2]
        pattern = "%Y%j%H%M%S"
        return datetime.strptime(filename, pattern)
