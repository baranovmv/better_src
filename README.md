# Sample Rate Converter that Rocs

These design choices lead to certain compromises:

* Time stamps in fixed point format can't hold rational values perfectly. For example, if we upsample
44.1 kHz to 96 kHz, it means the upsampling coefficient becomes 960/441, or the timestamp that holds
a position of an outgoing sample increments by 441/960, which doesn't fit perfectly in Q12.20 (it becomes
440.999450684 / 960). These tiny errors makes the signal to be sampled in a bit jittery time grid, which converts
to a phase modulation which could be represented as a raised noise floor:


  
Though it looks scary the absolute numbers are so low, it can't be perceived by a human being.  