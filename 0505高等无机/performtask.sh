#!/bin/bash

g09 ir2cl8-stable.gjf &
wait
g09 ir2cl8-45.gjf &
wait
g09 os2cl8-45.gjf &
wait
g09 os2cl8-stable.gjf &
wait
g09 ir2cl8-opt.gjf &
wait
g09 os2cl8-opt.gjf &
