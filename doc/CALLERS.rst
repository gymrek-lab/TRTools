Supported TR Callers
====================

TRTools currently supports 5 tandem repeat callers.
Here we introduce these callers and provide some basic specification of their functionality.
For more information on a caller, please see its website linked below.

+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
|                         |      AdVNTR_             | ExpansionHunter_        | GangSTR_                | HipSTR_                  | PopSTR_ (v2.0)          |
+=========================+==========================+=========================+=========================+==========================+=========================+
| Input Read Type         | Short Read or Long Read  | Short Read              | Short Read              | Short Read               | Short Read              |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
| Maximum Motif Size (bp) | 100                      | 6                       | 20                      | 6                        | 6                       |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
| Allele Size Limit       | Shorter than read length | Longer than read length | Longer than read length | Shorter than read length | Longer than read length |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
| Joint Calling           | No                       | No                      | No                      | Yes                      | Yes                     |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+
| Reference Based         | Yes                      | Yes                     | Yes                     | Yes                      | Yes                     |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------+-------------------------+

TRTools can be extended to support other callers that generate standard VCF files.
We welcome community contributions help support them. If that interests you, please 
see :ref:`Contributing` for more information.

..
    please ensure this list of links remains the same as the one in the main README

.. _AdVNTR: https://advntr.readthedocs.io/en/latest/
.. _ExpansionHunter: https://github.com/Illumina/ExpansionHunter
.. _GangSTR: https://github.com/gymreklab/gangstr
.. _HipSTR: https://hipstr-tool.github.io/HipSTR/
.. _PopSTR: https://github.com/DecodeGenetics/popSTR
 
