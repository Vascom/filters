#!/bin/bash

quartus_map --64bit --parallel=on --read_settings_files=on --write_settings_files=off leon3mp -c leon3mp &&
quartus_cdb --64bit --read_settings_files=on --write_settings_files=off leon3mp -c leon3mp --merge=on &&
quartus_fit --64bit --parallel=8 --read_settings_files=on --write_settings_files=off leon3mp -c leon3mp &&
quartus_sta --64bit --parallel=8 leon3mp -c leon3mp
quartus_asm --64bit --read_settings_files=on --write_settings_files=off leon3mp -c leon3mp
