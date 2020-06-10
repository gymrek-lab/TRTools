Supported TR Callers
====================

TR Tools currently supports 5 TR Callers. In this table, we introduce these callers and provide some basic specification of their functionality.

+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
|          Method         |          AdVNTR          | ExpansionHunter         | GangSTR                 | HipSTR                   | PopSTR2                 |
+=========================+==========================+=========================+=========================+==========================+=========================+
|  Input Sequencing Reads |  Short Read or Long Read | Short Read              | Short Read              | Short Read               | Short Read              |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
| Maximum Motif Size (bp) | 100                      | 6                       | 20                      | 6                        | 6                       |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
| Allele Size Limit       | Shorter than read length | Longer than read length | Longer than read length | Shorter than read length | Longer than read length |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
| Joint Calling           | No                       | No                      | No                      | Yes                      | Yes                     |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
| Reference Based         | Yes                      | Yes                     | Yes                     | Yes                      | Yes                     |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+

TR Tools can be extended to support other existing or new callers that generate a standard VCF files. We welcome contributions from the community to help expand the functionality of TR Tools.



