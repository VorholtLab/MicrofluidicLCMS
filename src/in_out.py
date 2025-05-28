# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 09:31:24 2023

@author: pkiefer
"""
import yaml
import os

script_dir = os.path.abspath(os.path.dirname(__file__))
data_folder = os.path.join(script_dir, '..', 'Data')


def load_config(path):
    with open(path, "r") as fp:
        data = fp.read()
    return yaml.safe_load(data)


