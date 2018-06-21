module mod_sym

  use mod_prec, only: rp, ip
  implicit none
  private

  real(kind=rp), allocatable, dimension(:,:), public :: xyzcom
  
  !
  ! symmetry params
  !
  ! TOLdist ........ Check of a distance.
  ! TOLeqvm ........ Comparison of matrices.
  ! TOLsng ......... Singularity of a 3x3 matrix.
  ! TOLisint ....... Closeness to an integer(kind=ip).
  ! TOLnull ........ Null angles and projections.
  ! TOLeigen ....... Eigenvalues of the orthogonal matrices.
  ! TOLdirty ....... Activate the algorithm changes to deal with
  !                  noisy data.
  ! inf_order ...... Infinity order will be transformed into this.
  ! TOLtriplet...... Used to select an initial triplet of atoms to
  !                  determine the symmetry matrices. This value
  !                  should be large (0.1 is the default). Increase
  !                  this value to ask for even better triplets.
  !
  logical, public :: TOLdirty
  integer(kind=ip), public :: inf_order
  real(kind=rp), public :: TOLdist, TOLeqvm, TOLsng, TOLisint, TOLnull
  real(kind=rp), public :: TOLeigen, TOLtriplet
  !
  ! The largest errors on the different tests.
  ! 
  ! ERRortho ....... Error on the orthogonality of the matrices.
  ! ERRsng ......... Best nonsingular triplet found.
  ! ERReigen ....... Error on the eigenvalues.
  !
  real(kind=rp), public :: ERRsng, ERReigen
  ! 
  ! nopsym ........ Number of symmetry operations.
  ! optype() ...... Symmetry operation types (see below).
  ! oporder() ..... Rotation order (the n in C_n^m or S_n^m).
  ! opm() ......... The m in C_n^m or S_n^m.
  ! opinv() ....... Inverse of an operation.
  ! opsym(i,,) .... The 3x3 matrix of the operation in the xyz rep.
  ! opaxis(i,) .... Rotation axis.
  ! opeuler(i,).... Euler angles for the rotation axis.
  ! opangle()...... Rotation angle (radians).
  ! opsymbol() .... Symbol for the operation.
  ! opproper() .... Proper/improper rotation (redundant).
  ! opisgener() ... Is this a generator for the group?
  ! optable(,)..... Multiplication (Cayley) table.
  ! linear_mol .... True if the molecule is linear.
  ! point_group ... Point group symbol.
  ! point_g2 ...... Point group symbol (short version).
  ! opmainord ..... Order of the main axis.
  ! 
  integer(kind=ip), parameter, public :: MOPSYM = 480
  logical, public  :: opproper(MOPSYM), opisgener(MOPSYM)
  integer(kind=ip), public  :: nopsym, optype(MOPSYM), oporder(MOPSYM)
  integer(kind=ip), public  :: opm(MOPSYM), opinv(MOPSYM)
  integer(kind=ip), public  :: optable(MOPSYM,MOPSYM), opmainord
  real(kind=rp), public :: opsym(MOPSYM,3,3)
  real(kind=rp), public :: opaxis(MOPSYM,3), opeuler(MOPSYM,3)
  real(kind=rp), public :: opangle(MOPSYM)
  character(len=10), public :: opsymbol(MOPSYM)
  character(len=20), public :: point_group
  character(len=8), public  :: point_g2
  !
  ! Matrices used to convert the molecule to a symmetry orientation:
  !
  real(kind=rp), public :: or_mat(3,3), or_imat(3,3)
  ! 
  ! Classes of symmetry operations:
  !
  integer(kind=ip), parameter :: MCLAS = 120, MINCLAS = 200
  integer(kind=ip), public :: nclas, opclas(MOPSYM)
  integer(kind=ip), public :: ninclas(MCLAS), iinclas(MCLAS,MINCLAS)
  !
  ! Sym operator types:
  !
  integer(kind=ip), parameter, public :: opt_identity     = 10
  integer(kind=ip), parameter, public :: opt_rotation     = 11
  integer(kind=ip), parameter, public :: opt_inversion    = 12
  integer(kind=ip), parameter, public :: opt_sigma        = 13
  integer(kind=ip), parameter, public :: opt_imp_rotation = 14
  integer(kind=ip), parameter, public :: opt_unknown      = 15
  !
  ! xcm .. zcm ...... center of mass coordinates.
  ! norbit .......... number of orbits.
  ! natorb() ........ number of atoms in an orbit.
  ! iatorb(,) ....... list of atoms in an orbit.
  ! orbZ() .......... atomic number of the atoms in this orbit.
  ! orbdis() ........ distance to the center for all
  !                   atoms in this orbit.
  ! orbmol() ........ orbit to which an atom belongs.
  ! molradius ....... maximum distance from the center to an atom.
  !
  integer(kind=ip), parameter, public :: MORBIT = 2000, MATORB = 300, MMOL1 = 300000
  integer(kind=ip), public  :: norbit, natorb(MORBIT), iatorb(MORBIT,MATORB)
  integer(kind=ip), public  :: orbZ(MORBIT), orbmol(MMOL1)
  real(kind=rp), public :: orbdis(MORBIT), molradius, xcm, ycm, zcm

  ! public routines
  !public :: sym1d, sym2d, sym3d, symgetgroup, syminit, symheader
  !public :: symback, symcleancoord, symclosure, symdeter, symeigen 
  !public :: symeuler, symimagic, syminv, symmatfill, symmatprod 
  !public :: symnoise, symopadd, symorb, symprint, symproject, sympurify
  !public :: symreplicate, symrwcoord, symtransform, symorient, symopchk

  ! aux logical data
  logical, public :: mol_linear, mol_planar, linear_mol

end module mod_sym
