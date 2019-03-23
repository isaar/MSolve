using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;


namespace ISAAR.MSolve.MultiscaleAnalysis.Interfaces
{
    /// <summary>
    /// Indicates the primary stress and strain tensor that will be used in the used multiscale homogenization scheme for the micro to macro transitions
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public interface IScaleTransitions_v2
    {
        double[] MacroToMicroTransition(Node_v2 boundaryNode, double[] MacroScaleVariable);
        double[] MicroToMacroTransition(INode boundaryNode, double[] MicroScaleVariable);
        void ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(Node_v2 boundaryNode,
            double[] MacroScaleVariable, Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements);
        void ImposeAppropriateConstraintsPerBoundaryNode(Model_v2 model, Node_v2 boundaryNode);
        void ImposeAppropriateAndRigidBodyConstraintsPerBoundaryNode(Model_v2 model_v2, Node_v2 boundaryNode, Dictionary<Node_v2, IList<DOFType>> RigidBodyNodeConstraints); //TODO: enopoihsh twn duo duplicate
        int PrescribedDofsPerNode(); // TODO: pithanws epistrofh kai to poioi einai me input sugkekrimeno node
        int MacroscaleVariableDimension();
    }
}
