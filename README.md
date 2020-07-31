# Extended-GC
Part 1： EGCI_CODE
Two neuron network Cases:
(a) Fig 3: HH_SCI_MAIN_5_28401_Addnoise_all.m
(b) Fig 4: HH_SCI_MAIN_5_280_Addnoise_all.m
(c) Fig 5: HH_SCI_MAIN_371_Addnoise_all.m
(d) Fig 11: HH_SCI_MAIN_5_3766_Addnoise_all.m

Three neuron network Cases:
(a) HH_Three_280290371:
f1=0.280 f2=0.290 f3=0.371
X->Z; Y->Z

(b) HH_Three_280:
f1=0.280 f2=0.282 f3=0.290
X->Y,X->Z

(c) HH_Three_280290371_2:
f1=0.280 f2=0.290 f3=0.371
X->Y; Y->Z

(d) HH_Three_280290371_Q:
Neutype=0: 表示 兴奋性和抑制性神经元的结构是 1 1 1, 结果在HH_Three_280290371文件夹
Neutype=1: 表示 兴奋性和抑制性神经元的结构是 1 1 0
Neutype=2: 表示 兴奋性和抑制性神经元的结构是 1 0 0

Neutype: 1 1 0
f1=0.280 f2=0.290 f3=0.371
X->Z; Y->Z

f1=0.280 f2=0.290 f3=0.371
NeuronType: 1  1  0
X->Y; Y->Z

Discussion Part
Embeding Dimension Compute
1) EGCI_Chen_EMD.m % Chen's paper's example
2) EGCI_Li_EMD.m  % Li's paper 's example


Supplementary Code
1) Three neuron network
1: causality_paper_GC_Three_HH_demo.m
2: causality_paper_GC_Three_HH_demo_xyycode.m

2) Two neuron network
1: causality_paper_GC_Two_HH_demo.m
2: causality_paper_GC_Two_HH_demo_xyycode.m


Part 2：H-H model simulation 
Through this program, voltage time series can obtained by evolving the HH neuronal network.


Part 3：Different_System_Final_Wolf
Wolf’s algorithm and compare with the LLE of the original HH neuronal network
using the standard algorithm
