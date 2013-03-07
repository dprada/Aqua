Kinetic Analysis
****************

1D trajectories
===============

bla, bla, bla

Ejemplos de committed:

El primer y ultimo indice solo pueden aparecer una vez.

aa=kinetic_1D_analysis([1,0,1,0,2])
aa.first_committed_passage_time(states=[1,2],commitment=[True,True])
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,True,True])
(array([ 2.]), array([ 1.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,False,True])
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,True,True])
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,False,True])
(array([ 2.,  4.]), array([ 1.,  1.]))


aa=kinetic_1D_analysis([1,1,1,2,1,1,1,0,1,2,1,1,3,2])
aa.first_committed_passage_time(states=[1,2],commitment=[True,True])
(array([ 1.,  2.,  3.]), array([ 2.,  1.,  1.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,True,True])
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,False,True])
(array([ 1.,  2.,  3.]), array([ 2.,  2.,  2.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,True,True])
(array([ 2.,  3.]), array([ 1.,  1.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,False,True])
(array([ 1.,  2.,  3.,  4.,  5.]), array([ 2.,  1.,  2.,  1.,  1.]))
aa.first_committed_passage_time(states=[1,4,2],commitment=[True,True,True])
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,4,2],commitment=[True,False,True])
(array([ 1.,  2.,  3.,  4.,  5.]), array([ 2.,  2.,  3.,  1.,  1.]))

aa=kinetic_1D_analysis([1,3,1,2,1,0,1,0,1,2,1,0,0,2])
aa.first_committed_passage_time(states=[1,2],commitment=[True,True])
(array([ 1.]), array([ 2.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,True,True])
(array([ 3.]), array([ 1.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,False,True])
(array([ 1.,  3.]), array([ 2.,  1.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,True,True])
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,False,True])
(array([ 1.,  3.,  5.]), array([ 2.,  2.,  1.]))
aa.first_committed_passage_time(states=[1,4,2],commitment=[True,True,True])
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,4,2],commitment=[True,False,True])
(array([ 1.,  3.,  5.]), array([ 2.,  3.,  1.]))


--- noreturn

aa=kinetic_1D_analysis([1,0,1,0,2])
aa.first_committed_passage_time(states=[1,2],commitment=[True,True],no_return=True)
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,True,True],no_return=True)
(array([ 2.]), array([ 1.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,False,True],no_return=True)
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,True,True],no_return=True)
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,False,True],no_return=True)
(array([ 2.]), array([ 1.]))


aa=kinetic_1D_analysis([1,1,1,2,1,1,1,0,1,2,1,1,3,2])
aa.first_committed_passage_time(states=[1,2],commitment=[True,True],no_return=True)
(array([ 1.,  2.,  3.]), array([ 2.,  1.,  1.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,True,True],no_return=True)
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,False,True],no_return=True)
(array([ 1.,  2.,  3.]), array([ 2.,  2.,  2.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,True,True],no_return=True)
(array([ 2.,  3.]), array([ 1.,  1.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,False,True],no_return=True)
(array([ 1.,  2.,  3.]), array([ 2.,  1.,  1.]))
aa.first_committed_passage_time(states=[1,4,2],commitment=[True,True,True],no_return=True)
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,4,2],commitment=[True,False,True],no_return=True)
(array([ 1.,  2.,  3.]), array([ 2.,  2.,  2.]))


aa=kinetic_1D_analysis([1,3,1,2,1,0,1,0,1,2,1,0,0,2])
aa.first_committed_passage_time(states=[1,2],commitment=[True,True],no_return=True)
(array([ 1.]), array([ 2.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,True,True],no_return=True)
(array([ 3.]), array([ 1.]))
aa.first_committed_passage_time(states=[1,0,2],commitment=[True,False,True],no_return=True)
(array([ 1.]), array([ 2.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,True,True],no_return=True)
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,3,2],commitment=[True,False,True],no_return=True)
(array([ 1.,  3.]), array([ 2.,  1.]))
aa.first_committed_passage_time(states=[1,4,2],commitment=[True,True,True],no_return=True)
(array([ 0.]), array([ 0.]))
aa.first_committed_passage_time(states=[1,4,2],commitment=[True,False,True],no_return=True)
(array([ 1.,  3.]), array([ 2.,  1.]))

