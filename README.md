# ContinuumRobotExamples

Continuum robots have elastic links which are capable of large-scale continuous deformations. This repo has example scripts to simulate continuum robots of various design paradigms. Each example is accompanied by a short write-up in PDF format.

![Blender rendering of a tendon-actuated robot](https://github.com/JohnDTill/ContinuumRobotExamples/raw/master/01_Statics/05_Tendon_Robot/LaTeX/fig/TendonRobotRender.jpg "Tendon-Driven Robot")

## Motivation

* Introduce continuum robot modeling
* Demonstrate advanced concepts such as dynamic simulation
* Disseminate published academic results

Continuum and soft robotics problems can have a steep learning curve. The mechanics of robot motions are described by differential equations. Controlling a robot may depend on writing a simulation loop to numerically solve for robot behavior at interactive rates. Although there is fundamental complexity which is difficult to generalize in software, this repo aims to provide an easier path with examples as learning aids.

The examples begin simple and build to present progressively more complicated models. Code is C++ with some examples also implemented in MATLAB. Qt is recommended to use with the C++ examples as the project organization files are in ".pro" Qt format, Qt's plotting capabilities are used to visualize models, and the event system is used to demonstrate interactive simulation.

![Interactive control of a parallel continuum robot](https://github.com/JohnDTill/ContinuumRobotExamples/raw/master/CSG_Control.gif "Interactive control of a parallel continuum robot")

## Prerequisites

* Qt (https://www.qt.io/download)
* PDF Viewer
* MATLAB (optional)
* Blender 2.7.9 (https://download.blender.org/release/Blender2.79/) (optional for creating procedural renderings. Note this is an **old version** due to breaking changes in scripting syntax.)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Contributing

Feel free to submit pull requests and use the issue tracker to start a discussion about any bugs you encounter. Please provide a description of your compiler and operating system for any software related bugs.

## Acknowledgements

Thanks to my Ph.D. advisor Caleb Rucker for his excellent work leading the [REACH lab](https://sites.google.com/site/danielcrucker/), and
thanks to the continuum robotics research community for creating many exciting robot designs and publishing great papers showing how one can use mathematical modeling to gain deeper understanding of physical systems.
