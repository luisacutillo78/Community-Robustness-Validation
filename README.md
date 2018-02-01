# Community-Robustness-Validation

In this repository we provide the basic R code and few examples to replicate results described in:

https://www.sciencedirect.com/science/article/pii/S0167947317302347

We provide an effective method to validate the robustness of a given community structure in a Network. 

Our main idea is summarized in the following abstract.

**Abstract** <br>
The large amount of work on community detection and its applications leaves unaddressed one important question: the statistical validation of the results. A methodology is presented that is able to clearly detect if the community structure found by some algorithms is statistically significant or is a result of chance, merely due to edge positions in the network. Given a community detection method and a network of interest, the proposal examines the stability of the partition recovered against random perturbations of the original graph structure. To address this issue, a perturbation strategy and a null model graph, which matches the original in some of its structural properties, but is otherwise a random graph, is specified. A set of procedures is built based on a special measure of clustering distance, namely Variation of Information, using tools set up for functional data analysis. The procedures determine whether the obtained clustering departs significantly from the null model. This strongly supports the robustness against perturbation of the algorithm used to identify the community structure. Results obtained with the proposed technique on simulated and real datasets are shown and discussed.

Please feel free to try it and provide us with a feedback!
