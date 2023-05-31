# VideoSelectionDEOUDAK

# Video Selection with Differential Evolution

We have 10 video files, each given with size and the play time.

![image](https://user-images.githubusercontent.com/22652266/235850052-ae9a0ab5-67d1-4cb2-bc76-9ad01d5fe156.png)

Need to store them in a DVD, to maximize the playing time while the total size of the files not excluding 4500 MB.

## Representation Method 

Binary Representation

## Fitness Function 

Fitness Function = (Σ weights*playtimes) * f
Where {f = 1, (Σ weights*cdsizes) <= 4500
f = 0.5, (Σ weights*cdsizes) > 4500}

## Mutation operator – Differential mutation

Here new mutant vector which is called donor vector is generated by adding weighted difference between 2 two vectors to another vector from population. These are distinct vectors and not the target vector. But since binary representation used here, the values generated with the above cannot be used directly (since it generates values other than 1 and 0). So as an assumption, values which are neither 1 nor 0 are assigned 1 or 0 randomly.

## Crossover operator – Uniform Crossover

Here uniform crossover is performed between donor vector generated from above mutation and target vector and generates a vector called trial vector.

## Parameters

N = 10
F = 1
CR = 0.1
Iterations = 100 (1 run)

## Information Obtained for 20 Runs

1st Run
bestChromoFitness = 609
ind = 1
Total size in mb = 4450
ans = 1   1   1   1   1   0   1   0   0   0

2nd Run
bestChromoFitness = 613
ind = 3
Total size in mb = 4475
ans = 0   1   0   0   0   1   1   1   1   0

3rd Run
bestChromoFitness = 614
ind = 10
Total size in mb = 4425
ans = 1   0   0   1   0   0   1   1   1   0

4th Run
bestChromoFitness = 603
ind = 1
Total size in mb = 4425
ans = 0   0   1   0   0   1   1   1   1   0

5th Run
bestChromoFitness = 589
ind = 4
Total size in mb = 4300
ans = 1   0   0   0   1   1   1   0   1   0

6th Run
bestChromoFitness = 614
ind = 8
Total size in mb = 4425
ans = 1   0   0   1   0   0   1   1   1   0
   
7th Run
bestChromoFitness = 614
ind = 9
ans = 4425
ans = 1   0   0   1   0   0   1   1   1   0
   
8th Run
bestChromoFitness = 613
ind = 3
Total size in mb = 4475
ans = 0   1   0   0   0   1   1   1   1   0

9th Run
bestChromoFitness = 614
ind = 8
Total size in mb = 4425
ans = 1   0   0   1   0   0   1   1   1   0
   
10th Run
bestChromoFitness = 614
ind = 1
Total size in mb = 4425
ans = 1   0   0   1   0   0   1   1   1   0
   
11th Run
bestChromoFitness = 613
ind = 1
Total size in mb = 4475
ans = 0   1   0   0   0   1   1   1   1   0
   
12th Run
bestChromoFitness = 611
ind = 4
Total size in mb = 4475
ans = 0   1   1   1   1   1   0   1   0   0
   
13th Run
bestChromoFitness = 614
ind = 2
Total size in mb = 4425
ans = 1   0   0   1   0   0   1   1   1   0

14th Run
bestChromoFitness = 614
ind = 1
Total size in mb = 4425
ans = 1   0   0   1   0   0   1   1   1   0
   
15th Run
bestChromoFitness = 614
ind = 1
Total size in mb = 4425
ans = 1   0   0   1   0   0   1   1   1   0
   
16th Run
bestChromoFitness = 613
ind = 5
Total size in mb = 4475
ans = 0   1   0   0   0   1   1   1   1   0

17th Run
bestChromoFitness = 603
ind = 2
Total size in mb = 4425
ans = 0   0   1   0   0   1   1   1   1   0
18th Run
bestChromoFitness = 614
ind = 6
Total size in mb = 4425
ans = 1   0   0   1   0   0   1   1   1   0

19th Run
bestChromoFitness = 607
ind = 4
Total size in mb= 4375
ans = 1   1   1   1   1   0   0   1   0   0
   
20th Run
bestChromoFitness = 614
ind = 2
Total size in mb = 4425
ans = 1   0   0   1   0   0   1   1   1   0

