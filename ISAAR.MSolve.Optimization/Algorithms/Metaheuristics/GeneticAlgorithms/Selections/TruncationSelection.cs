using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Optimization.Commons;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections
{
    public class TruncationSelection<T> : ISelectionStrategy<T>
    {
        private readonly double elitePercentage;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="elitePercentage">Legal values: (0, 1]. Recommended values: [0.1, 0.5]</param>
        public TruncationSelection(double elitePercentage = 0.5)
        {
            if (elitePercentage <= 0 || elitePercentage > 1)
            {
                throw new ArgumentException("The percentage of individuals that are eligible for selection must" 
                    + "belong to the interval (0, 1], but was: " + elitePercentage);
            }
            this.elitePercentage = elitePercentage;
        }

        public Individual<T>[][] Apply(Individual<T>[] population, int parentGroupsCount, int parentsPerGroup)
        {
            // Gather the parents that will be selected 
            int elitesCount = (int)Math.Round(elitePercentage * population.Length);
            if (elitesCount == 0) elitesCount = 1;
            Array.Sort(population);
            var selectionPool = new Individual<T>[elitesCount];
            Array.Copy(population, selectionPool, elitesCount);

            // Figure out how many copies of each selected parent will be used
            Individual<T>[] allParents = SelectEnoughElites(selectionPool, parentGroupsCount * parentsPerGroup);

            // Shuffle them to avoid grouping together identical parents or the fittest ones.
            RandomNumberGenerationUtilities.Shuffle<Individual<T>>(allParents);

            // Place the parents in groups
            int counter = 0;
            var parentGroups = new Individual<T>[parentGroupsCount][];
            for (int group = 0; group < parentGroupsCount; ++group)
            {
                parentGroups[group] = new Individual<T>[parentsPerGroup];
                for (int parent = 0; parent < parentsPerGroup; ++parent)
                {
                    parentGroups[group][parent] = allParents[counter];
                    ++counter;
                }
            }
            return parentGroups;
        }

        private Individual<T>[] SelectEnoughElites(Individual<T>[] selectedParents, int totalParentsCount)
        {
            // Sequentially select individuals starting from the fittest ones. Repeat while necessary.
            var allParents = new Individual<T>[totalParentsCount];
            int totalCounter = 0;
            int selectedCounter = 0;
            while (totalCounter < totalParentsCount)
            {
                allParents[totalCounter++] = selectedParents[selectedCounter++];
                if (selectedCounter == selectedParents.Length) selectedCounter = 0;
            }
            return allParents;
        }
    }
}
