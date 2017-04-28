
# DEM General Options
Dimension                        = 3
PeriodicDomainOption             = "OFF"
BoundingBoxOption                = "OFF"
BoundingBoxEnlargementFactor     = 1.0
AutomaticBoundingBoxOption       = "OFF"
BoundingBoxMaxX                  =  1.00000e+01
BoundingBoxMaxY                  =  1.00000e+01
BoundingBoxMaxZ                  =  1.00000e+01
BoundingBoxMinX                  = -1.00000e+01
BoundingBoxMinY                  = -1.00000e+01
BoundingBoxMinZ                  = -1.00000e+01

dem_inlet_option                 = 0
GravityX                         = 0.0
GravityY                         = 0.0
GravityZ                         = -9.81

EnergyCalculationOption          = 0
VelocityTrapOption               = 0
RotationOption                   = "ON"
CleanIndentationsOption          = "OFF"
DeltaOption                      = "Absolute"
SearchTolerance                  = 0.0001
AmplifiedSearchRadiusExtension   = 1e-4
MaxAmplificationRatioOfSearchRadius = 100
ModelDataInfo                    = "OFF"
VirtualMassCoefficient           = 1.0
RollingFrictionOption            = "OFF"
ComputeStressTensorOption        = "OFF"
PoissonEffectOption              = "ON"
ShearStrainParallelToBondOption  = "ON"
DontSearchUntilFailure           = "OFF"
ContactMeshOption                = "OFF"
OutputFileType                   = "Binary"
Multifile                        = "multiple_files"
ElementType                      = "SphericContPartDEMElement3D"

# Solution Strategy
IntegrationScheme                = "Symplectic_Euler"
AutomaticTimestep                = "OFF"
DeltaTimeSafetyFactor            = 1.0
MaxTimeStep                      = 5.0e-5
FinalTime                        = 1.0
ControlTime                      = 4.0
NeighbourSearchFrequency         = 10

# PostProcess Results
GraphExportFreq                  = 1e-3
VelTrapGraphExportFreq           = 1e0
OutputTimeStep                   = 1e-1
PostBoundingBox                  = 0
PostDisplacement                 = 1
PostVelocity                     = 1
# DEM only Results
PostTotalForces                  = 1
PostRigidElementForces           = 0
PostSkinSphere                   = 0
PostPoissonRatio                 = 0
PostRadius                       = 0
PostAngularVelocity              = 0
PostParticleMoment               = 0
PostEulerAngles                  = 0
PostRollingResistanceMoment      = 0
# FEM only Results
PostElasticForces                = 0
PostContactForces                = 0
PostTangentialElasticForces      = 0
PostShearStress                  = 0
PostPressure                     = 0
# FEM_wear only Results
PostNonDimensionalVolumeWear     = 0
PostNodalArea                    = 0
# Results on bond elements
PostStressStrainOption           = 0
PostContactSigma                 = 0
PostContactTau                   = 0
PostLocalContactForce            = 0
PostFailureCriterionState        = 0
PostContactFailureId             = 0
PostMeanContactArea              = 0
# Under revision
PostRHS                          = 0
PostDampForces                   = 0
PostAppliedForces                = 0
PostGroupId                      = 0
PostExportId                     = 0

#
problem_name="fibers_and_balls"
