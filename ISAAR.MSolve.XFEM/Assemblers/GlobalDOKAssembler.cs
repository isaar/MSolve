using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

/// <summary>
/// The matrix that will be "inverted" is a unified DOK matrix, where enriched dofs are numbered after all 
/// standard dofs. 
/// TODO: The enriched dof columns will have huge heights. A more sophisticated solver and matrix assembler are
/// needed. Also the global constrained submatrix must be sparse.
/// </summary>
namespace ISAAR.MSolve.XFEM.Assemblers
{
    class GlobalDOKAssembler: GlobalAssemblerBase<DOKSymmetricColMajor>
    {
        protected override DOKSymmetricColMajor InitializeGlobalUncontrainedMatrix(Model2D model, IDOFEnumerator dofEnumerator)
        {
            //TODO: empty or identity?
            return DOKSymmetricColMajor.CreateEmpty(dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount);
        }
    }
}
