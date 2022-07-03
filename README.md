Public Key Encryption with Conjunctive Keyword Search over NTRU Lattices
===========

This software is a proof-of-concept implementation of a PECKS scheme over NTRU lattices. The paper hasn't been published yet.

Warning
=======
This code is not to be considered secure, efficient or fully portable. Its purpose is not to be used for actual encryption, but to provide the research community a tool to verify, analyze and reproduce the statements made in our paper.

How to use?
===========

To modify the parameters, edit the values N0 and q0 in params.h.

To run on an Unix machine with g++:
```
$ make
$ ./PECKS
```

If GMP and NTL are not in a standard directory, you have to modify the CCFLAGS and LDFLAGS in the Makefile to indicate where they are.

