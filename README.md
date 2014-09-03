TTD
===

Tracking Training Day

## Setup

1. Setup a proper CMSSW release area and source its environments
```
cmsrel CMSSW_7_2_0_pre5
cd CMSSW_7_2_0_pre5
eval `scram r -sh`
```

1. To facilitate the access to your own github repository w/o
constantly typing your passphrase, do
```
eval `ssh-agent`
ssh-add
```
and digit your passphrase at the prompt

1. Setup an empty git repository for the CMSSW area you just created
```
git cms-init
```

1. Add code for the Tracking Training Day, 2nd edition
```
git submodule add git@github.com:/rovere/TTD TTD
```
