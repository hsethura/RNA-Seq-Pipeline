#!/bin/sh


mkdir gs_info

cd gs_info

cp /broad/IDP-Dx_work/nirmalya/gs_extract/client_secret.json .

/broad/IDP-Dx_work/nirmalya/pipeline/beta/quickstart.py -s 1fg8Ob825A4T9FfFZpmMogIaUl8MX4RdP5zbGM9RpGmE -t "Sample Information" -p Orna_New  --noauth_local_webserver

