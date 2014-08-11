#include <OpenSim/OpenSim.h>
#include "yaml-cpp/yaml.h"

int main()
{
    try {

        YAML::Node parameters = YAML::LoadFile("../../data/example_constants.yml");
        std::cout << parameters["mc"].as<double>() << std::endl;

        OpenSim::Model osimModel = OpenSim::Model();

        osimModel.setUseVisualizer(true);

        osimModel.setName("Gait2D");

        osimModel.setGravity(SimTK::Vec3(0, -parameters["g"].as<double>(), 0));

        OpenSim::Body &ground = osimModel.getGroundBody();

        /* Floor */
        /* The floor has a longitudinal degree of freedom and can be perturbed
         * by a longitudinal force or have a prescribed velocity. I'll treat
         * this as an infinite floor sliding under the user (just like a
         * treadmill). */

        double floorMass = 100;
        SimTK::Vec3 comLocInBody(0.0, 0.0, 0.0);
        SimTK::Inertia bodyInertia(10.0, 10.0, 1.0, 0.0, 0.0, 0.0);

        OpenSim::Body* floor = new OpenSim::Body("Floor",
                                                 floorMass,
                                                 comLocInBody,
                                                 bodyInertia);

        floor->addDisplayGeometry("box.vtp");
        floor->updDisplayer()->setScaleFactors(SimTK::Vec3(1.0, 0.05, 1.0));

        osimModel.addBody(floor);

        SimTK::Vec3 locationInParent(0.0, 0.0, 0.0);
        SimTK::Vec3 orientationInParent(0.0, 0.0, 0.0);
        SimTK::Vec3 locationInChild(0.0, 0.0, 0.0);
        SimTK::Vec3 orientationInChild(0.0, 0.0, 0.0);

        /* A SliderJoint translates along the x axis */
        OpenSim::SliderJoint *floorToGround =
            new OpenSim::SliderJoint("FloorToGround",
                                     ground,
                                     locationInParent,
                                     orientationInParent,
                                     *floor,
                                     locationInChild,
                                     orientationInChild,
                                     false);

        OpenSim::CoordinateSet &floorJoints = floorToGround->upd_CoordinateSet();

        floorJoints[0].setName("floor_translate_x");
        double translationRangeFloor[2] = {-0.1, 2.0 * 8.0 * 60.0};  // 2 m/s * 8 min * 60 sec/min
        floorJoints[0].setRange(translationRangeFloor);
        floorJoints[0].setDefaultValue(0.0);

        osimModel.addJoint(floorToGround);

        /* Trunk */
        /* The trunk (hip/torso/head/arms) can translate and rotate relative to
         * the fixed ground */

        double trunkLength = 2.0 *  parameters["ya"].as<double>();

        comLocInBody = SimTK::Vec3(parameters["xa"].as<double>(),
                                   parameters["ya"].as<double>(),
                                   0.0);

        bodyInertia = SimTK::Inertia(10.0, 10.0, parameters["ia"].as<double>(),
                                     0.0, 0.0, 0.0);

        OpenSim::Body* trunk = new OpenSim::Body("Trunk",
                                                 parameters["ia"].as<double>(),
                                                 comLocInBody,
                                                 bodyInertia);

        locationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::FreeJoint *trunkToFloor =
            new OpenSim::FreeJoint("TrunkToGround",
                                   ground,
                                   locationInParent,
                                   orientationInParent,
                                   *trunk,
                                   locationInChild,
                                   orientationInChild,
                                   false);

        /* This trunk translates in the x and y directions and rotates about z. */

        OpenSim::CoordinateSet &trunkJointCoords = trunkToFloor->upd_CoordinateSet();

        /* Rotation */

        trunkJointCoords[0].setName("qa_rotx");
        double xRotRangeTrunk[2] = {-SimTK::Pi, -SimTK::Pi};
        trunkJointCoords[0].setRange(xRotRangeTrunk);
        trunkJointCoords[0].setDefaultValue(SimTK::convertDegreesToRadians(0.0));
        trunkJointCoords[0].setDefaultLocked(true);

        trunkJointCoords[1].setName("qa_roty");
        double yRotRangeTrunk[2] = {-SimTK::Pi, SimTK::Pi};
        trunkJointCoords[1].setRange(yRotRangeTrunk);
        trunkJointCoords[1].setDefaultValue(SimTK::convertDegreesToRadians(0.0));
        trunkJointCoords[1].setDefaultLocked(true);

        trunkJointCoords[2].setName("qa");
        double zRotRangeTrunk[2] = {-SimTK::Pi, SimTK::Pi};
        trunkJointCoords[2].setRange(zRotRangeTrunk);
        trunkJointCoords[2].setDefaultValue(SimTK::convertDegreesToRadians(0.0));

        /* Translation */

        trunkJointCoords[3].setName("qax");
        double xTranRangeTrunk[2] = {-10.0, 100.0};
        trunkJointCoords[3].setRange(xTranRangeTrunk);
        trunkJointCoords[3].setDefaultValue(0.0);

        trunkJointCoords[4].setName("qay");
        double yTranRangeTrunk[2] = {-1.0, 2.0};
        trunkJointCoords[4].setRange(yTranRangeTrunk);
        trunkJointCoords[4].setDefaultValue(1.0); // TODO : This should be some nominal height, probably from measured data.

        trunkJointCoords[5].setName("qaz");
        double zTranRangeTrunk[2] = {-1.0, 1.0};
        trunkJointCoords[5].setRange(zTranRangeTrunk);
        trunkJointCoords[5].setDefaultValue(0.0);
        trunkJointCoords[5].setDefaultLocked(true);

        /* Visualization Geometry */

        trunk->addDisplayGeometry("sphere.vtp");
        trunk->updDisplayer()->setScaleFactors(SimTK::Vec3(trunkLength / 10.0, trunkLength, trunkLength / 10.0));
        trunk->updDisplayer()->translate(SimTK::Vec3(0.0, trunkLength / 2.0, 0.0));

        osimModel.addBody(trunk);
        osimModel.addJoint(trunkToFloor);

        /* Right Thigh: Body B */

        comLocInBody = SimTK::Vec3(parameters["xb"].as<double>(),
                                   parameters["yb"].as<double>(),
                                   0.0);

        bodyInertia = SimTK::Inertia(10.0, 10.0, parameters["ib"].as<double>(),
                                     0.0, 0.0, 0.0);

        OpenSim::Body* rightThigh = new OpenSim::Body("RightThigh",
                                                      parameters["mb"].as<double>(),
                                                      comLocInBody,
                                                      bodyInertia);

        double rightThighLength = parameters["lb"].as<double>();

        locationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);

        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::PinJoint *rightThighToTrunk =
            new OpenSim::PinJoint("RightThighToTrunk",
                                  *trunk,
                                  locationInParent,
                                  orientationInParent,
                                  *rightThigh,
                                  locationInChild,
                                  orientationInChild,
                                  false);

        OpenSim::CoordinateSet &rightThighJoints = rightThighToTrunk->upd_CoordinateSet();

        rightThighJoints[0].setName("qb");
        double zRotRangeRightThigh[2] = {-SimTK::convertDegreesToRadians(100.0),
                                         SimTK::convertDegreesToRadians(100.0)};
        rightThighJoints[0].setRange(zRotRangeRightThigh);
        rightThighJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(30.0));

        rightThigh->addDisplayGeometry("sphere.vtp");
        rightThigh->updDisplayer()->setScaleFactors(SimTK::Vec3(rightThighLength / 10.0,
                                                                rightThighLength,
                                                                rightThighLength / 10.0));
        rightThigh->updDisplayer()->translate(SimTK::Vec3(0.0, -rightThighLength / 2.0, 0.0));

        osimModel.addBody(rightThigh);
        osimModel.addJoint(rightThighToTrunk);

        /* Right Shank : Body C */

        comLocInBody = SimTK::Vec3(parameters["xc"].as<double>(),
                                   parameters["yc"].as<double>(),
                                   0.0);
        bodyInertia = SimTK::Inertia(10.0, 10.0, parameters["ic"].as<double>(),
                                     0.0, 0.0, 0.0);

        OpenSim::Body* rightShank = new OpenSim::Body("Right Shank",
                                                      parameters["mc"].as<double>(),
                                                      comLocInBody,
                                                      bodyInertia);

        double rightShankLength = parameters["lc"].as<double>();

        locationInParent = SimTK::Vec3(0.0, -rightThighLength, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);

        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::PinJoint *rightShankToThigh =
            new OpenSim::PinJoint("RightShankToThigh",
                                  *rightThigh,
                                  locationInParent,
                                  orientationInParent,
                                  *rightShank,
                                  locationInChild,
                                  orientationInChild,
                                  false);

        OpenSim::CoordinateSet &rightShankJoints = rightShankToThigh->upd_CoordinateSet();

        rightShankJoints[0].setName("qc");
        double zRotRangeRightShank[2] = {-SimTK::convertDegreesToRadians(100.0), 0.0};
        rightShankJoints[0].setRange(zRotRangeRightShank);
        rightShankJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(-30.0));

        rightShank->addDisplayGeometry("sphere.vtp");
        rightShank->updDisplayer()->setScaleFactors(SimTK::Vec3(rightShankLength / 10.0,
                                                                rightShankLength,
                                                                rightShankLength / 10.0));
        rightShank->updDisplayer()->translate(SimTK::Vec3(0.0, -rightShankLength / 2.0, 0.0));

        osimModel.addBody(rightShank);
        osimModel.addJoint(rightShankToThigh);

        /* Right Foot : Body D */

        comLocInBody = SimTK::Vec3(parameters["xd"].as<double>(),
                                   parameters["yd"].as<double>(),
                                   0.0);

        bodyInertia = SimTK::Inertia(10.0, 10.0, parameters["id"].as<double>(),
                                     0.0, 0.0, 0.0);

        OpenSim::Body* rightFoot = new OpenSim::Body("Right Foot",
                                                      parameters["md"].as<double>(),
                                                      comLocInBody,
                                                      bodyInertia);

        double rightToeDistance = parameters["txd"].as<double>();
        double rightHeelDistance = parameters["hxd"].as<double>();
        double rightFootDepth = parameters["fyd"].as<double>();

        double footLength = rightToeDistance - rightHeelDistance;

        locationInParent = SimTK::Vec3(0.0, -rightShankLength, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);

        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::PinJoint *rightFootToShank =
            new OpenSim::PinJoint("RightFootToShank",
                                  *rightShank,
                                  locationInParent,
                                  orientationInParent,
                                  *rightFoot,
                                  locationInChild,
                                  orientationInChild,
                                  false);

        OpenSim::CoordinateSet &rightFootJoints = rightFootToShank->upd_CoordinateSet();

        rightFootJoints[0].setName("qd");
        double zRotRangeRightFoot[2] = {-SimTK::convertDegreesToRadians(100.0), 0.0};
        rightFootJoints[0].setRange(zRotRangeRightFoot);
        rightFootJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(-30.0));

        rightFoot->addDisplayGeometry("sphere.vtp");
        rightFoot->updDisplayer()->setScaleFactors(SimTK::Vec3(footLength,
                                                               footLength / 10.0,
                                                               footLength / 10.0));

        osimModel.addBody(rightFoot);
        osimModel.addJoint(rightFootToShank);

        /* Left Thigh: Body E */

        comLocInBody = SimTK::Vec3(parameters["xe"].as<double>(),
                                   parameters["ye"].as<double>(),
                                   0.0);

        bodyInertia = SimTK::Inertia(10.0, 10.0, parameters["ie"].as<double>(),
                                     0.0, 0.0, 0.0);

        OpenSim::Body* leftThigh = new OpenSim::Body("LeftThigh",
                                                     parameters["me"].as<double>(),
                                                     comLocInBody,
                                                     bodyInertia);

        double leftThighLength = parameters["le"].as<double>();

        locationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);

        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::PinJoint *leftThighToTrunk =
            new OpenSim::PinJoint("LeftThighToTrunk",
                                  *trunk,
                                  locationInParent,
                                  orientationInParent,
                                  *leftThigh,
                                  locationInChild,
                                  orientationInChild,
                                  false);

        OpenSim::CoordinateSet &leftThighJoints = leftThighToTrunk->upd_CoordinateSet();

        leftThighJoints[0].setName("qe");
        double zRotRangeLeftThigh[2] = {-SimTK::convertDegreesToRadians(100.0),
                                         SimTK::convertDegreesToRadians(100.0)};
        leftThighJoints[0].setRange(zRotRangeLeftThigh);
        leftThighJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(-10.0));

        leftThigh->addDisplayGeometry("sphere.vtp");
        leftThigh->updDisplayer()->setScaleFactors(SimTK::Vec3(leftThighLength / 10.0,
                                                               leftThighLength,
                                                               leftThighLength / 10.0));
        leftThigh->updDisplayer()->translate(SimTK::Vec3(0.0, -leftThighLength / 2.0, 0.0));

        osimModel.addBody(leftThigh);
        osimModel.addJoint(leftThighToTrunk);

        /* Left Shank : Body F */

        comLocInBody = SimTK::Vec3(parameters["xf"].as<double>(),
                                   parameters["yf"].as<double>(),
                                   0.0);
        bodyInertia = SimTK::Inertia(10.0, 10.0, parameters["if"].as<double>(),
                                     0.0, 0.0, 0.0);

        OpenSim::Body* leftShank = new OpenSim::Body("LeftShank",
                                                     parameters["mf"].as<double>(),
                                                     comLocInBody,
                                                     bodyInertia);

        double leftShankLength = parameters["lf"].as<double>();

        locationInParent = SimTK::Vec3(0.0, -leftThighLength, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);

        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::PinJoint *leftShankToThigh =
            new OpenSim::PinJoint("LeftShankToThigh",
                                  *leftThigh,
                                  locationInParent,
                                  orientationInParent,
                                  *leftShank,
                                  locationInChild,
                                  orientationInChild,
                                  false);

        OpenSim::CoordinateSet &leftShankJoints = leftShankToThigh->upd_CoordinateSet();

        leftShankJoints[0].setName("qf");
        double zRotRangeLeftShank[2] = {-SimTK::convertDegreesToRadians(100.0), 0.0};
        leftShankJoints[0].setRange(zRotRangeLeftShank);
        leftShankJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(-30.0));

        leftShank->addDisplayGeometry("sphere.vtp");
        leftShank->updDisplayer()->setScaleFactors(SimTK::Vec3(leftShankLength / 10.0,
                                                               leftShankLength,
                                                               leftShankLength / 10.0));
        leftShank->updDisplayer()->translate(SimTK::Vec3(0.0, -leftShankLength / 2.0, 0.0));

        osimModel.addBody(leftShank);
        osimModel.addJoint(leftShankToThigh);

        /* Left Foot : Body G */

        comLocInBody = SimTK::Vec3(parameters["xg"].as<double>(),
                                   parameters["yg"].as<double>(),
                                   0.0);

        bodyInertia = SimTK::Inertia(10.0, 10.0, parameters["ig"].as<double>(),
                                     0.0, 0.0, 0.0);

        OpenSim::Body* leftFoot = new OpenSim::Body("Left Foot",
                                                    parameters["mg"].as<double>(),
                                                    comLocInBody,
                                                    bodyInertia);

        double leftToeDistance = parameters["txg"].as<double>();
        double leftHeelDistance = parameters["hxg"].as<double>();
        double leftFootDepth = parameters["fyg"].as<double>();

        double leftFootLength = leftToeDistance - leftHeelDistance;

        locationInParent = SimTK::Vec3(0.0, -leftShankLength, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);

        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::PinJoint *leftFootToShank =
            new OpenSim::PinJoint("LeftFootToShank",
                                  *leftShank,
                                  locationInParent,
                                  orientationInParent,
                                  *leftFoot,
                                  locationInChild,
                                  orientationInChild,
                                  false);

        OpenSim::CoordinateSet &leftFootJoints = leftFootToShank->upd_CoordinateSet();

        leftFootJoints[0].setName("qg");
        double zRotRangeLeftFoot[2] = {-SimTK::convertDegreesToRadians(100.0), 0.0};
        leftFootJoints[0].setRange(zRotRangeLeftFoot);
        leftFootJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(-30.0));

        leftFoot->addDisplayGeometry("sphere.vtp");
        leftFoot->updDisplayer()->setScaleFactors(SimTK::Vec3(leftFootLength,
                                                              leftFootLength / 10.0,
                                                              leftFootLength / 10.0));

        osimModel.addBody(leftFoot);
        osimModel.addJoint(leftFootToShank);

        /* Contact Geometry */

        OpenSim::ContactHalfSpace* floorContact =
            new OpenSim::ContactHalfSpace(SimTK::Vec3(0.0, 0.0, 0.0),
                                          SimTK::Vec3(0.0, 0.0, -SimTK::Pi / 2.0),
                                          *floor,
                                          "FloorContact");
        osimModel.addContactGeometry(floorContact);

        double contactSphereRadius = 0.05;

        /* Hip */

        SimTK::Vec3 rightHipLocationInTrunk(0.0, 0.0, 0.0);

        OpenSim::ContactSphere* rightHipContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       rightHipLocationInTrunk,
                                       *trunk,
                                       "RightHipContact");

        osimModel.addContactGeometry(rightHipContact);

        SimTK::Vec3 leftHipLocationInTrunk(0.0, 0.0, 0.0);

        OpenSim::ContactSphere* leftHipContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       leftHipLocationInTrunk,
                                       *trunk,
                                       "LeftHipContact");
        osimModel.addContactGeometry(leftHipContact);

        /* Knee */

        SimTK::Vec3 rightKneeLocationInThigh(0.0, 0.0, 0.0);

        OpenSim::ContactSphere* rightKneeContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       rightKneeLocationInThigh,
                                       *rightShank,
                                       "RightKneeContact");

        osimModel.addContactGeometry(rightKneeContact);

        SimTK::Vec3 leftKneeLocationInThigh(0.0, 0.0, 0.0);

        OpenSim::ContactSphere* leftKneeContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       leftKneeLocationInThigh,
                                       *leftShank,
                                       "LeftKneeContact");

        osimModel.addContactGeometry(leftKneeContact);

        /* Ankle */

        SimTK::Vec3 rightAnkleLocationInThigh(0.0, 0.0, 0.0);

        OpenSim::ContactSphere* rightAnkleContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       rightAnkleLocationInThigh,
                                       *rightShank,
                                       "RightAnkleContact");

        osimModel.addContactGeometry(rightAnkleContact);

        SimTK::Vec3 leftAnkleLocationInThigh(0.0, 0.0, 0.0);

        OpenSim::ContactSphere* leftAnkleContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       leftAnkleLocationInThigh,
                                       *leftShank,
                                       "LeftAnkleContact");

        osimModel.addContactGeometry(leftAnkleContact);

        /* Feet */
        /* A sphere is offset from the heel and toe points. This is slightly
         * different than the pygait2d model because it uses a point
         * definition instead of a sphere. */

        double heelContactSphereRadius = 0.06;
        double toeContactSphereRadius = 0.04;

        SimTK::Vec3 rightHeelLocationInFoot(rightHeelDistance + heelContactSphereRadius,
                                            rightFootDepth + heelContactSphereRadius,
                                            0.0);

        OpenSim::ContactSphere* rightHeelContact =
            new OpenSim::ContactSphere(heelContactSphereRadius,
                                       rightHeelLocationInFoot,
                                       *rightFoot,
                                       "RightHeelContact");

        osimModel.addContactGeometry(rightHeelContact);

        SimTK::Vec3 rightToeLocationInFoot(rightToeDistance - toeContactSphereRadius,
                                           rightFootDepth + toeContactSphereRadius,
                                           0.0);

        OpenSim::ContactSphere* rightToeContact =
            new OpenSim::ContactSphere(toeContactSphereRadius,
                                       rightToeLocationInFoot,
                                       *rightFoot,
                                       "RightToeContact");

        osimModel.addContactGeometry(rightToeContact);

        SimTK::Vec3 leftHeelLocationInFoot(leftHeelDistance + heelContactSphereRadius,
                                           leftFootDepth + heelContactSphereRadius,
                                            0.0);

        OpenSim::ContactSphere* leftHeelContact =
            new OpenSim::ContactSphere(heelContactSphereRadius,
                                       leftHeelLocationInFoot,
                                       *leftFoot,
                                       "LeftHeelContact");

        osimModel.addContactGeometry(leftHeelContact);

        SimTK::Vec3 leftToeLocationInFoot(leftToeDistance - toeContactSphereRadius,
                                          leftFootDepth + toeContactSphereRadius,
                                           0.0);

        OpenSim::ContactSphere* leftToeContact =
            new OpenSim::ContactSphere(toeContactSphereRadius,
                                       leftToeLocationInFoot,
                                       *leftFoot,
                                       "LeftToeContact");

        osimModel.addContactGeometry(leftToeContact);


        /* Contact Forces */

        // TODO : I may not be mapping these correctly as I couldn't find the
        // equation used to compute this force anywhere.

        double stiffness = parameters["kc"].as<double>();
        double dissipation = parameters["cc"].as<double>();
        double staticFriction = parameters["mu"].as<double>();
        double dynamicFriction = parameters["mu"].as<double>();
        double viscosity = parameters["vs"].as<double>();

        /* Hip */
        OpenSim::HuntCrossleyForce::ContactParameters *rightHipContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        rightHipContactParams->addGeometry("RightHipContact");
        rightHipContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* rightHipForce =
            new OpenSim::HuntCrossleyForce(rightHipContactParams);
        rightHipForce->setName("RightHipForce");

        osimModel.addForce(rightHipForce);

        OpenSim::HuntCrossleyForce::ContactParameters *leftHipContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        leftHipContactParams->addGeometry("LeftHipContact");
        leftHipContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* leftHipForce =
            new OpenSim::HuntCrossleyForce(leftHipContactParams);
        leftHipForce->setName("LeftHipForce");

        osimModel.addForce(leftHipForce);

        /* Knee */
        OpenSim::HuntCrossleyForce::ContactParameters *rightKneeContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        rightKneeContactParams->addGeometry("RightKneeContact");
        rightKneeContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* rightKneeForce =
            new OpenSim::HuntCrossleyForce(rightKneeContactParams);
        rightKneeForce->setName("RightKneeForce");

        osimModel.addForce(rightKneeForce);

        OpenSim::HuntCrossleyForce::ContactParameters *leftKneeContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        leftKneeContactParams->addGeometry("LeftKneeContact");
        leftKneeContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* leftKneeForce =
            new OpenSim::HuntCrossleyForce(leftKneeContactParams);
        leftKneeForce->setName("LeftKneeForce");

        osimModel.addForce(leftKneeForce);

        /* Ankle */

        OpenSim::HuntCrossleyForce::ContactParameters *rightAnkleContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);

        rightAnkleContactParams->addGeometry("RightAnkleContact");
        rightAnkleContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* rightAnkleForce =
            new OpenSim::HuntCrossleyForce(rightAnkleContactParams);
        rightAnkleForce->setName("RightAnkleForce");

        osimModel.addForce(rightAnkleForce);

        OpenSim::HuntCrossleyForce::ContactParameters *leftAnkleContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);

        leftAnkleContactParams->addGeometry("LeftAnkleContact");
        leftAnkleContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* leftAnkleForce =
            new OpenSim::HuntCrossleyForce(leftAnkleContactParams);
        leftAnkleForce->setName("LeftAnkleForce");

        osimModel.addForce(leftAnkleForce);

        /* Foot */

        /* right heel */

        OpenSim::HuntCrossleyForce::ContactParameters *rightHeelContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        rightHeelContactParams->addGeometry("RightHeelContact");
        rightHeelContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* rightHeelForce =
            new OpenSim::HuntCrossleyForce(rightHeelContactParams);
        rightHeelForce->setName("RightHeelForce");

        osimModel.addForce(rightHeelForce);

        /* right toe */

        OpenSim::HuntCrossleyForce::ContactParameters *rightToeContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        rightToeContactParams->addGeometry("RightToeContact");
        rightToeContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* rightToeForce =
            new OpenSim::HuntCrossleyForce(rightToeContactParams);
        rightToeForce->setName("RightToeForce");

        osimModel.addForce(rightToeForce);

        /* left heel */

        OpenSim::HuntCrossleyForce::ContactParameters *leftHeelContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        leftHeelContactParams->addGeometry("LeftHeelContact");
        leftHeelContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* leftHeelForce =
            new OpenSim::HuntCrossleyForce(leftHeelContactParams);
        leftHeelForce->setName("LeftHeelForce");

        osimModel.addForce(leftHeelForce);

        /* left toe */

        OpenSim::HuntCrossleyForce::ContactParameters *leftToeContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        leftToeContactParams->addGeometry("LeftToeContact");
        leftToeContactParams->addGeometry("FloorContact");

        OpenSim::HuntCrossleyForce* leftToeForce =
            new OpenSim::HuntCrossleyForce(leftToeContactParams);
        leftToeForce->setName("LeftToeForce");

        osimModel.addForce(leftToeForce);

        /* Coordinate Limits */
        /* The pygait2d and algait2d models do not yet have coordinate limit
         * forces. Plus I want the identification to find these, not impose
         * them, so they are commented out. */

        /*

        double limitStiffness = 1e6;
        double limitDamping = 1e5;
        double limitTransition = 5.0; // degrees

        */

        /* Hip */
        /*
        OpenSim::CoordinateLimitForce* rightHipLimit =
            new OpenSim::CoordinateLimitForce("RHip_rz",
                                              zRotRangeRightThigh[0],
                                              limitStiffness,
                                              zRotRangeRightThigh[1],
                                              limitStiffness,
                                              limitDamping,
                                              limitTransition,
                                              false);
        osimModel.addForce(rightHipLimit);

        OpenSim::CoordinateLimitForce* leftHipLimit =
            new OpenSim::CoordinateLimitForce("LHip_rz",
                                              zRotRangeLeftThigh[0],
                                              limitStiffness,
                                              zRotRangeLeftThigh[1],
                                              limitStiffness,
                                              limitDamping,
                                              limitTransition,
                                              false);
        osimModel.addForce(leftHipLimit);
        */

        /* Knee */
        /*
        OpenSim::CoordinateLimitForce* rightKneeLimit =
            new OpenSim::CoordinateLimitForce("RKnee_rz",
                                              zRotRangeRightShank[0],
                                              limitStiffness,
                                              zRotRangeRightShank[1],
                                              limitStiffness,
                                              limitDamping,
                                              limitTransition,
                                              false);
        osimModel.addForce(rightKneeLimit);

        OpenSim::CoordinateLimitForce* leftKneeLimit =
            new OpenSim::CoordinateLimitForce("LKnee_rz",
                                              zRotRangeLeftShank[0],
                                              limitStiffness,
                                              zRotRangeLeftShank[1],
                                              limitStiffness,
                                              limitDamping,
                                              limitTransition,
                                              false);
        osimModel.addForce(leftKneeLimit);

        */

        /* Generalized Loads */
        /* This adds joint torques for the ankle, knee, and hip. */

        double maxTorque = 1000.0;

        OpenSim::CoordinateActuator* rightHipTorque = new OpenSim::CoordinateActuator("qb");
        rightHipTorque->setName("RightHipTorque");
        rightHipTorque->setMinControl(-maxTorque);
        rightHipTorque->setMaxControl(maxTorque);
        osimModel.addForce(rightHipTorque);

        OpenSim::CoordinateActuator* rightKneeTorque = new OpenSim::CoordinateActuator("qc");
        rightKneeTorque->setName("RightKneeTorque");
        rightKneeTorque->setMinControl(-maxTorque);
        rightKneeTorque->setMaxControl(maxTorque);
        osimModel.addForce(rightKneeTorque);

        OpenSim::CoordinateActuator* rightAnkleTorque = new OpenSim::CoordinateActuator("qd");
        rightAnkleTorque->setName("RightAnkleTorque");
        rightAnkleTorque->setMinControl(-maxTorque);
        rightAnkleTorque->setMaxControl(maxTorque);
        osimModel.addForce(rightAnkleTorque);

        OpenSim::CoordinateActuator* leftHipTorque = new OpenSim::CoordinateActuator("qe");
        leftHipTorque->setName("LeftHipTorque");
        leftHipTorque->setMinControl(-maxTorque);
        leftHipTorque->setMaxControl(maxTorque);
        osimModel.addForce(leftHipTorque);

        OpenSim::CoordinateActuator* leftKneeTorque = new OpenSim::CoordinateActuator("qf");
        leftKneeTorque->setName("LeftKneeTorque");
        leftKneeTorque->setMinControl(-maxTorque);
        leftKneeTorque->setMaxControl(maxTorque);
        osimModel.addForce(leftKneeTorque);

        OpenSim::CoordinateActuator* leftAnkleTorque = new OpenSim::CoordinateActuator("qg");
        leftAnkleTorque->setName("LeftAnkleTorque");
        leftAnkleTorque->setMinControl(-maxTorque);
        leftAnkleTorque->setMaxControl(maxTorque);
        osimModel.addForce(leftAnkleTorque);

        /* Prescribed Motion for the ground */
        /* I'm going to move the ground under the walker */
        // SimTK::PrescribedMotion  this is a constraint

        /* Serialize Model */

        osimModel.print("Gait2D.osim");

        /* Simulate */


        // Configure the model.
        SimTK::State& state = osimModel.initSystem();

        // Add display geometry.
        osimModel.updMatterSubsystem().setShowDefaultGeometry(true);
        SimTK::Visualizer& viz = osimModel.updVisualizer().updSimbodyVisualizer();
        viz.setBackgroundColor(SimTK::Vec3(1, 1, 1));

        // Simulate.
        SimTK::RungeKuttaMersonIntegrator integrator(osimModel.getSystem());
        OpenSim::Manager manager(osimModel, integrator);
        manager.setInitialTime(0);
        manager.setFinalTime(10.0);
        manager.integrate(state);

    }
    catch (OpenSim::Exception ex)
    {
        std::cout << ex.getMessage() << std::endl;
        return 1;
    }
    catch (SimTK::Exception::Base ex)
    {
        std::cout << ex.getMessage() << std::endl;
        return 1;
    }
    catch (std::exception ex)
    {
        std::cout << ex.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cout << "UNRECOGNIZED EXCEPTION" << std::endl;
    }
    std::cout << "OpenSim example completed successfully" << std::endl;
    std::cout << "Press return to continue" << std::endl;
    std::cin.get();
    return 0;
}
