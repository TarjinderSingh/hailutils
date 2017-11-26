#!/usr/bin/python

import os
import sys
import re
import logging
from hail import *
from pprint import pprint

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)