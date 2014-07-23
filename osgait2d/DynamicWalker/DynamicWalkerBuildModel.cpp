#include <OpenSim/OpenSim.h>
#include "yaml-cpp/yaml.h"

int main()
{
    try {
        YAML::Node parameters = YAML::LoadFile("/home/moorepants/src/gait2d/data/example_constants.yml");
        std::cout << parameters["mc"].as<double>() << std::endl;

        double pelvisWidth = 0.20;
        double thighLength = 0.40;
        double shankLength = 0.435;
        double mass = 1;

        OpenSim::Model osimModel = OpenSim::Model();

        osimModel.setUseVisualizer(true);

        osimModel.setName("DynamicWalkerModel");

        osimModel.setGravity(SimTK::Vec3(0, -9.81, 0));
        OpenSim::Body &ground = osimModel.getGroundBody();

        /* Platform */

        SimTK::Vec3 comLocInBody(0.0, 0.0, 0.0);
        SimTK::Inertia bodyInertia(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

        OpenSim::Body* platform = new OpenSim::Body("Platform",
                                                    mass,
                                                    comLocInBody,
                                                    bodyInertia);

        platform->addDisplayGeometry("box.vtp");
        platform->updDisplayer()->setScaleFactors(SimTK::Vec3(1.0, 0.05, 1.0));

        osimModel.addBody(platform);

        SimTK::Vec3 locationInParent(0.0, 0.0, 0.0);
        SimTK::Vec3 orientationInParent(0.0, 0.0, 0.0);
        SimTK::Vec3 locationInChild(0.0, 0.0, 0.0);
        SimTK::Vec3 orientationInChild(0.0, 0.0, 0.0);

        OpenSim::PinJoint *platformToGround =
            new OpenSim::PinJoint("PlatformToGround",
                                  ground,
                                  locationInParent,
                                  orientationInParent,
                                  *platform,
                                  locationInChild,
                                  orientationInChild,
                                  false);

        OpenSim::CoordinateSet &platformJoints = platformToGround->upd_CoordinateSet();

        platformJoints[0].setName("platform_rz");
        double rotRangePlatform[2] = {-SimTK::Pi / 2.0, 0.0};
        platformJoints[0].setRange(rotRangePlatform);
        platformJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(-10.0));
        platformJoints[0].setDefaultLocked(true);

        osimModel.addJoint(platformToGround);

        /* Pelvis */

        comLocInBody = SimTK::Vec3(0.0, 0.0, 0.0);
        bodyInertia = SimTK::Inertia(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

        OpenSim::Body* pelvis = new OpenSim::Body("Pelvis", mass, comLocInBody, bodyInertia);

        locationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::FreeJoint *pelvisToPlatform =
            new OpenSim::FreeJoint("PelvisToPlatform",
                                   *platform,
                                   locationInParent,
                                   orientationInParent,
                                   *pelvis,
                                   locationInChild,
                                   orientationInChild,
                                   false);

        OpenSim::CoordinateSet &pelvisJointCoords = pelvisToPlatform->upd_CoordinateSet();


        /* This pelvis only translates in the x and y directions. */

        pelvisJointCoords[0].setName("pelvis_rx");
        double xRotRangePelvis[2] = {-SimTK::Pi, -SimTK::Pi};
        pelvisJointCoords[0].setRange(xRotRangePelvis);
        pelvisJointCoords[0].setDefaultValue(SimTK::convertDegreesToRadians(0.0));
        pelvisJointCoords[0].setDefaultLocked(true);

        pelvisJointCoords[1].setName("pelvis_ry");
        double yRotRangePelvis[2] = {-SimTK::Pi, SimTK::Pi};
        pelvisJointCoords[1].setRange(yRotRangePelvis);
        pelvisJointCoords[1].setDefaultValue(SimTK::convertDegreesToRadians(0.0));
        pelvisJointCoords[1].setDefaultLocked(true);

        pelvisJointCoords[2].setName("pelvis_rz");
        double zRotRangePelvis[2] = {-SimTK::Pi, SimTK::Pi};
        pelvisJointCoords[2].setRange(zRotRangePelvis);
        pelvisJointCoords[2].setDefaultValue(SimTK::convertDegreesToRadians(0.0));
        pelvisJointCoords[2].setDefaultLocked(true);

        pelvisJointCoords[3].setName("pelvis_tx");
        double xTranRangePelvis[2] = {-10.0, 10.0};
        pelvisJointCoords[3].setRange(xTranRangePelvis);
        pelvisJointCoords[3].setDefaultValue(0.0);

        pelvisJointCoords[4].setName("pelvis_ty");
        double yTranRangePelvis[2] = {-1.0, 2.0};
        pelvisJointCoords[4].setRange(yTranRangePelvis);
        pelvisJointCoords[4].setDefaultValue(1.0);

        pelvisJointCoords[5].setName("pelvis_tz");
        double zTranRangePelvis[2] = {-1.0, 1.0};
        pelvisJointCoords[5].setRange(zTranRangePelvis);
        pelvisJointCoords[5].setDefaultValue(0.0);
        pelvisJointCoords[5].setDefaultLocked(true);

        pelvis->addDisplayGeometry("sphere.vtp");
        pelvis->updDisplayer()->setScaleFactors(SimTK::Vec3(pelvisWidth / 2.0, pelvisWidth / 2.0, pelvisWidth));

        osimModel.addBody(pelvis);
        osimModel.addJoint(pelvisToPlatform);

        /* Left Thigh */

        comLocInBody = SimTK::Vec3(0.0, 0.0, 0.0);
        bodyInertia = SimTK::Inertia(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

        OpenSim::Body* leftThigh = new OpenSim::Body("LeftThigh", mass,
                                                     comLocInBody, bodyInertia);

        locationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::PinJoint *leftThighToPelvis =
            new OpenSim::PinJoint("LeftThighToPelvis",
                                  *pelvis,
                                  locationInParent,
                                  orientationInParent,
                                  *leftThigh,
                                  locationInChild,
                                  orientationInChild,
                                  false);

        OpenSim::CoordinateSet &leftThighJoints = leftThighToPelvis->upd_CoordinateSet();

        leftThighJoints[0].setName("LHip_rz");
        double zRotRangeLeftThigh[2] = {-SimTK::convertDegreesToRadians(100.0), SimTK::convertDegreesToRadians(100.0)};
        leftThighJoints[0].setRange(zRotRangeLeftThigh);
        leftThighJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(-10.0));

        leftThigh->addDisplayGeometry("sphere.vtp");
        leftThigh->updDisplayer()->setScaleFactors(SimTK::Vec3(thighLength / 10.0, thighLength, thighLength / 10.0));

        osimModel.addBody(leftThigh);
        osimModel.addJoint(leftThighToPelvis);

        /* Right Thigh */

        comLocInBody = SimTK::Vec3(0.0, 0.0, 0.0);
        bodyInertia = SimTK::Inertia(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

        OpenSim::Body* rightThigh = new OpenSim::Body("RightThigh", mass,
                                                      comLocInBody, bodyInertia);

        locationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
        locationInChild = SimTK::Vec3(0.0, 0.0, 0.0);
        orientationInChild = SimTK::Vec3(0.0, 0.0, 0.0);

        OpenSim::PinJoint *rightThighToPelvis =
            new OpenSim::PinJoint("RightThighToPelvis",
                                  *pelvis,
                                  locationInParent,
                                  orientationInParent,
                                  *rightThigh,
                                  locationInChild,
                                  orientationInChild,
                                  false);

        OpenSim::CoordinateSet &rightThighJoints = rightThighToPelvis->upd_CoordinateSet();

        rightThighJoints[0].setName("RHip_rz");
        double zRotRangeRightThigh[2] = {-SimTK::convertDegreesToRadians(100.0), SimTK::convertDegreesToRadians(100.0)};
        rightThighJoints[0].setRange(zRotRangeRightThigh);
        rightThighJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(30.0));

        rightThigh->addDisplayGeometry("sphere.vtp");
        rightThigh->updDisplayer()->setScaleFactors(SimTK::Vec3(thighLength / 10.0, thighLength, thighLength / 10.0));

        osimModel.addBody(rightThigh);
        osimModel.addJoint(rightThighToPelvis);

        /* Left Shank */

        comLocInBody = SimTK::Vec3(0.0, 0.0, 0.0);
        bodyInertia = SimTK::Inertia(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

        OpenSim::Body* leftShank = new OpenSim::Body("LeftShank", mass,
                                                     comLocInBody, bodyInertia);

        locationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
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

        leftShankJoints[0].setName("LKnee_rz");
        double zRotRangeLeftShank[2] = {-SimTK::convertDegreesToRadians(100.0), 0.0};
        leftShankJoints[0].setRange(zRotRangeLeftShank);
        leftShankJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(-30.0));

        leftShank->addDisplayGeometry("sphere.vtp");
        leftShank->updDisplayer()->setScaleFactors(SimTK::Vec3(shankLength / 10.0, shankLength, shankLength / 10.0));

        osimModel.addBody(leftShank);
        osimModel.addJoint(leftShankToThigh);

        /* Right Shank */

        comLocInBody = SimTK::Vec3(0.0, 0.0, 0.0);
        bodyInertia = SimTK::Inertia(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

        OpenSim::Body* rightShank = new OpenSim::Body("RightShank", mass,
                                                      comLocInBody, bodyInertia);

        locationInParent = SimTK::Vec3(0.0, 0.0, 0.0);
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

        rightShankJoints[0].setName("RKnee_rz");
        double zRotRangeRightShank[2] = {-SimTK::convertDegreesToRadians(100.0), 0.0};
        rightShankJoints[0].setRange(zRotRangeRightShank);
        rightShankJoints[0].setDefaultValue(SimTK::convertDegreesToRadians(-30.0));

        rightShank->addDisplayGeometry("sphere.vtp");
        rightShank->updDisplayer()->setScaleFactors(SimTK::Vec3(shankLength / 10.0, shankLength, shankLength / 10.0));

        osimModel.addBody(rightShank);
        osimModel.addJoint(rightShankToThigh);

        /* Contact Geometry */
        OpenSim::ContactHalfSpace* platformContact =
            new OpenSim::ContactHalfSpace(SimTK::Vec3(0.0, 0.0, 0.0),
                                          SimTK::Vec3(0.0, 0.0, -SimTK::Pi / 2.0),
                                          *platform,
                                          "PlatformContact");
        osimModel.addContactGeometry(platformContact);

        double contactSphereRadius = 0.05;

        /* Hip */

        SimTK::Vec3 rightHipLocationInPelvis(0.0, 0.0, pelvisWidth / 2.0);
        OpenSim::ContactSphere* rightHipContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       rightHipLocationInPelvis,
                                       *pelvis,
                                       "RHipContact");
        osimModel.addContactGeometry(rightHipContact);

        SimTK::Vec3 leftHipLocationInPelvis(0.0, 0.0, -pelvisWidth / 2.0);
        OpenSim::ContactSphere* leftHipContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       leftHipLocationInPelvis,
                                       *pelvis,
                                       "LHipContact");
        osimModel.addContactGeometry(leftHipContact);

        /* Knee */

        SimTK::Vec3 rightKneeLocationInThigh(0.0, thighLength / 2.0, pelvisWidth / 2.0);
        OpenSim::ContactSphere* rightKneeContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       rightKneeLocationInThigh,
                                       *rightThigh,
                                       "RKneeContact");
        osimModel.addContactGeometry(rightKneeContact);

        SimTK::Vec3 leftKneeLocationInThigh(0.0, thighLength / 2.0, -pelvisWidth / 2.0);
        OpenSim::ContactSphere* leftKneeContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       leftKneeLocationInThigh,
                                       *leftThigh,
                                       "LKneeContact");
        osimModel.addContactGeometry(leftKneeContact);

        /* Feet */

        SimTK::Vec3 rightFootLocationInThigh(0.0, shankLength / 2.0, pelvisWidth / 2.0);
        OpenSim::ContactSphere* rightFootContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       rightFootLocationInThigh,
                                       *rightShank,
                                       "RFootContact");
        osimModel.addContactGeometry(rightFootContact);

        SimTK::Vec3 leftFootLocationInThigh(0.0, shankLength / 2.0, -pelvisWidth / 2.0);
        OpenSim::ContactSphere* leftFootContact =
            new OpenSim::ContactSphere(contactSphereRadius,
                                       leftFootLocationInThigh,
                                       *leftShank,
                                       "LFootContact");
        osimModel.addContactGeometry(leftFootContact);

        /* Contact Forces */

        double stiffness = 1e7;
        double dissipation = 0.1;
        double staticFriction = 0.6;
        double dynamicFriction = 0.4;
        double viscosity = 0.01;

        /* Hip */
        OpenSim::HuntCrossleyForce::ContactParameters *rightHipContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        rightHipContactParams->addGeometry("RHipContact");
        rightHipContactParams->addGeometry("PlatformContact");

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
        leftHipContactParams->addGeometry("LHipContact");
        leftHipContactParams->addGeometry("PlatformContact");

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
        rightKneeContactParams->addGeometry("RKneeContact");
        rightKneeContactParams->addGeometry("PlatformContact");

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
        leftKneeContactParams->addGeometry("LKneeContact");
        leftKneeContactParams->addGeometry("PlatformContact");

        OpenSim::HuntCrossleyForce* leftKneeForce =
            new OpenSim::HuntCrossleyForce(leftKneeContactParams);
        leftKneeForce->setName("LeftKneeForce");

        osimModel.addForce(leftKneeForce);

        /* Foot */
        OpenSim::HuntCrossleyForce::ContactParameters *rightFootContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        rightFootContactParams->addGeometry("RFootContact");
        rightFootContactParams->addGeometry("PlatformContact");

        OpenSim::HuntCrossleyForce* rightFootForce =
            new OpenSim::HuntCrossleyForce(rightFootContactParams);
        rightFootForce->setName("RightFootForce");

        osimModel.addForce(rightFootForce);

        OpenSim::HuntCrossleyForce::ContactParameters *leftFootContactParams =
            new OpenSim::HuntCrossleyForce::ContactParameters(stiffness,
                                                              dissipation,
                                                              staticFriction,
                                                              dynamicFriction,
                                                              viscosity);
        leftFootContactParams->addGeometry("LFootContact");
        leftFootContactParams->addGeometry("PlatformContact");

        OpenSim::HuntCrossleyForce* leftFootForce =
            new OpenSim::HuntCrossleyForce(leftFootContactParams);
        leftFootForce->setName("LeftFootForce");

        osimModel.addForce(leftFootForce);

        /* Coordinate Limits */

        double limitStiffness = 1e6;
        double limitDamping = 1e5;
        double limitTransition = 5.0; // degrees

        /* Hip */
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

        /* Knee */
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

        /* Output */

        osimModel.print("DynamicWalkerModel.osim");

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
