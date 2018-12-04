using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: time logging must be refactored
namespace ISAAR.MSolve.Analyzers
{
    public class EquivalentLoadsAssembler_v2 : IEquivalentLoadsAssembler_v2
    {
        private ISubdomain_v2 subdomain;
        private IElementMatrixProvider elementProvider;

        public EquivalentLoadsAssembler_v2(ISubdomain_v2 subdomain, IElementMatrixProvider elementProvider)
        {
            this.subdomain = subdomain;
            this.elementProvider = elementProvider;
        }

        public IVector GetEquivalentNodalLoads(IVectorView solution, double constraintScalingFactor) 
        {
            //var times = new Dictionary<string, TimeSpan>();
            //var totalStart = DateTime.Now;
            //times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            //times.Add("element", TimeSpan.Zero);
            //times.Add("addition", TimeSpan.Zero);

            //var subdomainEquivalentNodalForces = new double[subdomain.Forces.Length];
            var subdomainEquivalentForces = Vector.CreateZero(subdomain.DofOrdering.NumFreeDofs);
            foreach (Element element in subdomain.Elements)
            {
                //var elStart = DateTime.Now;
                IMatrix elementK = elementProvider.Matrix(element).LegacyToNewMatrix();

                //double[] localSolution = subdomain.CalculateElementNodalDisplacements(element, solution);
                //double[] localdSolution = subdomain.CalculateElementIcrementalConstraintDisplacements(element, constraintScalingFactor);
                Vector localdSolution = subdomain.CalculateElementIncrementalConstraintDisplacements(element, constraintScalingFactor);

                //var equivalentNodalForces = new double[localSolution.Length];
                //ElementK.Multiply(new Vector(localdSolution), equivalentNodalForces);
                var elementEquivalentForces = elementK.Multiply(localdSolution);

                //subdomain.AddLocalVectorToGlobal(element, equivalentNodalForces, subdomainEquivalentNodalForces);
                subdomain.DofOrdering.AddVectorElementToSubdomain(element, elementEquivalentForces, subdomainEquivalentForces);

                //times["addition"] += DateTime.Now - elStart;
            }

            //var totalTime = DateTime.Now - totalStart;

            return subdomainEquivalentForces;
        }        
    }
}
