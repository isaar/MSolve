using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.Logging
{
    public class TotalLoadsDisplacementsPerIncrementLog_v1
    {
        private readonly Subdomain subdomain;
        //private readonly Dictionary<INode, HashSet<DOFType>> monitorDofs = new Dictionary<INode, HashSet<DOFType>>();
        private readonly INode monitorNode;
        private readonly DOFType monitorDof;
        private readonly string outputFile;

        public TotalLoadsDisplacementsPerIncrementLog_v1(Subdomain subdomain, int maxIncrementsExpected,
            INode monitorNode, DOFType monitorDof, string outputFile)
        {
            this.subdomain = subdomain;
            this.monitorNode = monitorNode;
            this.monitorDof = monitorDof;

            this.outputFile = outputFile;
        }

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

        public void LogTotalDataForIncrement(int incrementNumber, int currentIterationNumber, double errorNorm,
            IVector totalDisplacements, IVector totalInternalForces)
        {
            int subdomainDofIdx = subdomain.NodalDOFsDictionary[monitorNode.ID][monitorDof]; 

            // If all subdomains use the same file, then we need to open it in append mode. 
            //TODO: Also that will not work in parallel for many subdomains.
            using (var writer = new StreamWriter(outputFile, true)) // append mode to continue from previous increment
            {
                writer.Write($"{incrementNumber}, {currentIterationNumber}, {errorNorm}");
                writer.WriteLine($", {totalDisplacements[subdomainDofIdx]}, {totalInternalForces[subdomainDofIdx]}");
            }
        }
    }
}
