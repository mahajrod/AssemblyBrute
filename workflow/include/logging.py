#!/usr/bin/env python
import logging
import logging.config

from pathlib import Path

SPACE_PER_TAB = 4
TAB = " " * SPACE_PER_TAB

LOW_LVL_DEBUG_SCRIPT_NUM = 12
LOW_LVL_DEBUG_SCRIPT_NAME = "LOW_DBG_SCR"
LOW_LVL_DEBUG_SCRIPT_METHOD_NAME = LOW_LVL_DEBUG_SCRIPT_NAME.lower()

DEBUG_SCRIPT_NUM = 15
DEBUG_SCRIPT_NAME = "DBG_SCR"
DEBUG_SCRIPT_METHOD_NAME = DEBUG_SCRIPT_NAME.lower()

WARNING_SCRIPT_NUM = 25
WARNING_SCRIPT_NAME = "WARN_SCR"
WARNING_SCRIPT_METHOD_NAME = WARNING_SCRIPT_NAME.lower()

# Adding custom level to capture debug messages related to the script only

def logForLevelDBG_SCR(self, message, *args, **kwargs):
    if self.isEnabledFor(DEBUG_SCRIPT_NUM):
        self._log(DEBUG_SCRIPT_NUM, message, args, **kwargs)

def logToRootDBG_SCR(message, *args, **kwargs):
    logging.log(DEBUG_SCRIPT_NUM, message, *args, **kwargs)

def logForLevelLOW_LVL_DBG_SCR(self, message, *args, **kwargs):
    if self.isEnabledFor(LOW_LVL_DEBUG_SCRIPT_NUM):
        self._log(LOW_LVL_DEBUG_SCRIPT_NUM, message, args, **kwargs)

def logToRootLOW_LVL_DBG_SCR(message, *args, **kwargs):
    logging.log(LOW_LVL_DEBUG_SCRIPT_NUM, message, *args, **kwargs)

def logForLevelWARNING_SCR(self, message, *args, **kwargs):
    if self.isEnabledFor(WARNING_SCRIPT_NUM):
        self._log(WARNING_SCRIPT_NUM, message, args, **kwargs)

def logToRootWARNING_SCR(message, *args, **kwargs):
    logging.log(WARNING_SCRIPT_NUM, message, *args, **kwargs)

def filter_maker(level):
    level = getattr(logging, level)

    def filter(record):
        return record.levelno <= level

    return filter

FILTER_REGISTRY = {"warnings_and_below": filter_maker("WARNING")}

def get_filter(name):
    return FILTER_REGISTRY[name]


def config_logger(main_config):
    out_dir_path = Path(main_config["out_dir"])
    init_stage_log_prefix = main_config["init_stage_log_prefix"]


    logging.addLevelName(LOW_LVL_DEBUG_SCRIPT_NUM, LOW_LVL_DEBUG_SCRIPT_NAME)
    setattr(logging, LOW_LVL_DEBUG_SCRIPT_NAME, LOW_LVL_DEBUG_SCRIPT_NUM)
    setattr(logging.getLoggerClass(), LOW_LVL_DEBUG_SCRIPT_METHOD_NAME, logForLevelLOW_LVL_DBG_SCR)
    setattr(logging, LOW_LVL_DEBUG_SCRIPT_METHOD_NAME, logToRootLOW_LVL_DBG_SCR)

    logging.addLevelName(DEBUG_SCRIPT_NUM, DEBUG_SCRIPT_NAME)
    setattr(logging, DEBUG_SCRIPT_NAME, DEBUG_SCRIPT_NUM)
    setattr(logging.getLoggerClass(), DEBUG_SCRIPT_METHOD_NAME, logForLevelDBG_SCR)
    setattr(logging, DEBUG_SCRIPT_METHOD_NAME, logToRootDBG_SCR)

    logging.addLevelName(WARNING_SCRIPT_NUM, WARNING_SCRIPT_NAME)
    setattr(logging, WARNING_SCRIPT_NAME, WARNING_SCRIPT_NUM)
    setattr(logging.getLoggerClass(), WARNING_SCRIPT_METHOD_NAME, logForLevelWARNING_SCR)
    setattr(logging, WARNING_SCRIPT_METHOD_NAME, logToRootWARNING_SCR)

    log_config_dict = {
                       "version": 1,
                       "disable_existing_loggers": False,
                       "formatters": {
                                      "main": {
                                               "format": "%(asctime)s:%(levelname)-8s %(message)s",
                                               "datefmt": "%Y-%b-%d-%H:%M:%S",
                                              },
                                      "file": {
                                               "format": "%(asctime)s:%(levelname)-14s %(message)s",
                                               "datefmt": "%Y-%b-%d-%H:%M:%S",
                                              }
                                     },
                       "filters": {
                                   "warnings_and_below": {
                                                         "()" : get_filter,
                                                         "name": "warnings_and_below"
                                                         },
                                   },
                       "handlers": {
                           "stdout": {
                                      "class": "logging.StreamHandler",
                                      "level": "INFO",
                                      "formatter": "main",
                                      "stream": "ext://sys.stdout",
                                      "filters": ["warnings_and_below"]
                                      },
                           "stderr": {
                                     "class": "logging.StreamHandler",
                                     "level": "WARNING",
                                     "formatter": "main",
                                     "stream": "ext://sys.stderr"
                                      },
                           "debug_file": { # doesn't include debug messages from imported modules
                                         "class": "logging.FileHandler",
                                         "level": "DBG_SCR",
                                         "formatter": "file",
                                         "filename": out_dir_path / f"{init_stage_log_prefix}.debug_file.log",
                                         "mode": "w"
                                         },
                           "low_lvl_debug_file": { # doesn't include debug messages from imported modules
                                                 "class": "logging.FileHandler",
                                                 "level": "LOW_DBG_SCR",
                                                 "formatter": "file",
                                                 "filename": out_dir_path / f"{init_stage_log_prefix}.low_lvl_debug_file.log",
                                                 "mode": "w"
                                                 },
                           "full_file": { #includes debug messages from imported modules
                                         "class": "logging.FileHandler",
                                         "formatter": "file",
                                         "filename": out_dir_path / f"{init_stage_log_prefix}.total.log",
                                         "mode": "w"
                                        },
                       },
                       "loggers": {
                                   "brute_logger": {
                                                   "level": "DEBUG",
                                                   "handlers": [
                                                                "stderr",
                                                                "stdout",
                                                                "debug_file",
                                                                "low_lvl_debug_file",
                                                                "full_file"
                                                               ],
                                                   "propagate": False,
                                                   },

                                   },
                       }


    logging.config.dictConfig(log_config_dict)
    logger = logging.getLogger("brute_logger")

    return logger
