using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: This class should only extract data. How to output them (print in .txt, .xlsx, etc) should be done by different or
//      child classes
//TODO: Should this matrix write the results periodically (e.g. when a buffer fills up) instead of each increment?
//TODO: Extend it to extract data for many dofs simultaneously. It should also handles subdomains itself.
//TODO: If the analyzer ends abruptly (e.g. snap-through point in load control), that should be written here.
//TODO: Perhaps the file should not be opened and closed at each increment. Instead it should stay open, but then it should be 
//      disposed properly.
namespace ISAAR.MSolve.Logging
{
    public class TotalLoadsDisplacementsPerIncrementLog
    {
        private readonly ISubdomain_v2 subdomain;
        //private readonly Dictionary<INode, HashSet<DOFType>> monitorDofs = new Dictionary<INode, HashSet<DOFType>>();
        private readonly INode monitorNode;
        private readonly DOFType monitorDof;
        private readonly string outputFile;

        //TODO: These will be useful if I want to store them instead of immediately print them.
        //private readonly List<int> storedIncrementNumbers; //TODO: not sure if needed.
        //private readonly List<int> storedLastIterationNumbers;
        //private readonly List<double> storedErrorNorms;
        //private readonly List<double> storedTotalDisplacements;
        //private readonly List<double> storedTotalInternalForces;

        public TotalLoadsDisplacementsPerIncrementLog(ISubdomain_v2 subdomain, int maxIncrementsExpected, 
            INode monitorNode, DOFType monitorDof, string outputFile)
        {
            this.subdomain = subdomain;
            //this.storedIncrementNumbers = new List<int>(maxIncrementsExpected);
            //this.storedLastIterationNumbers = new List<int>(maxIncrementsExpected);
            //this.storedErrorNorms = new List<double>(maxIncrementsExpected);
            //this.storedTotalDisplacements = new List<double>(maxIncrementsExpected);
            //this.storedTotalInternalForces = new List<double>(maxIncrementsExpected);

            this.monitorNode = monitorNode;
            this.monitorDof = monitorDof;

            this.outputFile = outputFile;
        }

        /// <summary>
        /// Writes the header.
        /// </summary>
        public void Initialize()
        {
            // If all subdomains use the same file, then we need to open it in append mode. 
            //TODO: Also that will not work in parallel for many subdomains.
            using (var writer = new StreamWriter(outputFile, false)) // do not append, since this is a new analysis
            {
                // Header
                writer.Write("Increment, Iteration, ResidualNorm");
                writer.Write($", Total displacement (Node {monitorNode.ID} - dof {monitorDof})");
                writer.WriteLine($", Total internal force (Node {monitorNode.ID} - dof {monitorDof})");
            }
        }

        /// <summary>
        /// This also writes to the output file.
        /// </summary>
        /// <param name="totalDisplacements">
        /// The total displacements (start till current iteration of current increment) of the subdomain.
        /// </param>
        /// <param name="totalInternalForces">
        /// The total internal right hand side forces (start till current iteration of current increment) of the subdomain.
        /// </param>
        public void LogTotalDataForIncrement(int incrementNumber, int currentIterationNumber, double errorNorm,
            IVectorView totalDisplacements, IVectorView totalInternalForces)
        {
            int subdomainDofIdx = subdomain.DofOrdering.FreeDofs[monitorNode, monitorDof]; //TODO: SHould this be cached?

            // If all subdomains use the same file, then we need to open it in append mode. 
            //TODO: Also that will not work in parallel for many subdomains.
            using (var writer = new StreamWriter(outputFile, true)) // append mode to continue from previous increment
            {
                writer.Write($"{incrementNumber}, {currentIterationNumber}, {errorNorm}");
                writer.WriteLine($", {totalDisplacements[subdomainDofIdx]}, {totalInternalForces[subdomainDofIdx]}");
            }
        }

        //public void LogTotalDataForIncrement(int incrementNumber, int currentIterationNumber, double errorNorm,
        //    IVectorView totalDisplacements, IVectorView totalInternalForces)
        //{
        //    storedIncrementNumbers.Add(incrementNumber);
        //    storedLastIterationNumbers.Add(currentIterationNumber);
        //    storedErrorNorms.Add(errorNorm);

        //    int globalDofIdx = model.GlobalDofOrdering.GlobalFreeDofs[monitorNode, monitorDof];
        //    storedTotalDisplacements.Add(totalDisplacements[globalDofIdx]);
        //    storedTotalInternalForces.Add(totalInternalForces[globalDofIdx]);
        //}

        //public void Monitor(INode node, DOFType dof)
        //{
        //    bool nodeExists = monitorDofs.TryGetValue(node, out HashSet<DOFType> dofsOfNode);
        //    if (!nodeExists) dofsOfNode = new HashSet<DOFType>();
        //    dofsOfNode.Add(dof);
        //}
    }
}
