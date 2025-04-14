To generate tables showing the partial charges for the side chain analogues, the log files after training, along with the log file based on ESP-derived charges train_ff-ESP.log, are provided and will be used as input for the charges.py script.

- tune_ff-elec-p.log contains the training log for the electrostatic interaction energies using the point charge (PC) model.

- tune_ff-elec-g.log corresponds to the Gaussian charge (GC) model.

- tune_ff-elec-v.log covers the GC+PGV model (a non-polarizable Gaussian model with a virtual site).

- tune_ff-elec-h.log represents the PC+GVS model (a point charge model with a polarizable Gaussian-distributed shell).

The following log files correspond to models trained on both electrostatic and induction energies:

- tune_ff-allelec-i.log for the point charge (PC) model,

- tune_ff-allelec-s.log for the Gaussian charge (GC) model, and

- tune_ff-allelec-n.log for the GC+PGV model.
