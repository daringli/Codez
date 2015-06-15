12-06-2015 10:00
Codex_particle_model::Rsum should work now.
Excluded Ef=0, it now starts at Ef=dE, and that seems to have done the trick.

12-06-2015
Note: due to a rewritting Codex_particle_model::Rsum, the code currently calculates decay probabilities weirdly. It could be reverted back to an older (less optimized) version, but I plan this fix this the first thing I do next week.
