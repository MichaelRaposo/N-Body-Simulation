# N-Body Simulation of Celestial Bodies
## Multibody Jool Simulation
### Version 1 Notes

## Multibody Kerbol Simulation
Built upon Multibody Jool Simulation with emphasis on scalability and fixing form issues
### Version 1 Notes
- Moved Preallocation and Initilizations to beginning of document where possible
- Consistent Index, all arrays are same length with each index referring to same body
- Added Sphere Of Influence (SOI) to body property array
- Changed orbits array to match size of body_prop array
- added reference body to matrices in orbits array to indicate which body the orbital elements are in reference to (i.e. the parent body)
- Made various simulation constants seperate variables for ease of changes
- Cleaned up figure setup
- The trajectory and labels of moons (satelites) will not render until its distance from parent body is greater than SOI (i.e. it has left parent bodies SOI)
- Added progress tracker for longer simulations

### Version 2 Notes
- Changed simulation forloops to eliminate redundant calculations
- Greatly improved runtime

### Version 3 Notes


## Multibody Solar System Simulation
### Version 1 Notes
