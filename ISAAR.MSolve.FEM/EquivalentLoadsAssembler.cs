using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

//TODO: The subdomain should be passed as an argument in the required method, instead of being composed.
namespace ISAAR.MSolve.FEM
{
    public class EquivalentLoadsAssembler : IEquivalentLoadsAssembler
    {
        private Subdomain subdomain;
        private Discretization.Interfaces.IElementMatrixProvider elementProvider;

        public EquivalentLoadsAssembler(Subdomain subdomain, Discretization.Interfaces.IElementMatrixProvider elementProvider)
        {
            this.subdomain = subdomain;
            this.elementProvider = elementProvider;
        }

        public IVector GetEquivalentNodalLoads(IVector solution, double constraintScalingFactor) 
        {
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            var subdomainEquivalentNodalForces = new double[subdomain.Forces.Length];
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is Interfaces.IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);

                double[] localSolution = subdomain.CalculateElementNodalDisplacements(element, solution);
                double[] localdSolution = subdomain.CalculateElementIcrementalConstraintDisplacements(element, constraintScalingFactor);

                var equivalentNodalForces = new double[localSolution.Length];
                ElementK.Multiply(new Vector(localdSolution), equivalentNodalForces);
                subdomain.AddLocalVectorToGlobal(element, equivalentNodalForces, subdomainEquivalentNodalForces);

                times["addition"] += DateTime.Now - elStart;
            }

            var totalTime = DateTime.Now - totalStart;

            return new Vector(subdomainEquivalentNodalForces);
        }        
    }
}
