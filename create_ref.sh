# Linux Command Cheatsheet (Bash Script)

# General commands
# Select Interpreter (if running entire script)
#! /bin/bash

# View and change working directory
echo "Changed working directory from:"
pwd
cd /home/ncsteven/analyze/1_packages/HiVE
echo "to:"
pwd

# SQLite
## Create a new database in the current working directory
sqlite3 HiVE.db