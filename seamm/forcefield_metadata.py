metadata = {
    'charges':
        {
            'equation': ['I'],
            'constants': [
                ('Q', 'e'),
            ],
            'topology':
                {
                    'type': 'charges',
                    'n_atoms': 1,
                    'symmetry': 'none',
                    'fill': 0,
                    'flip': 0
                }
        },
    'bond_increments':
        {
            'equation': ['delta'],
            'constants': [
                ('deltaij', 'e'),
                ('deltaji', 'e'),
            ],
            'topology':
                {
                    'type': 'increment',
                    'n_atoms': 2,
                    'symmetry': 'like_bond',
                    'fill': 0,
                    'flip': 1
                }
        },
    'quadratic_bond':
        {
            'equation': ['K2*(R-R0)^2'],
            'constants': [
                ('R0', 'angstrom'),
                ('K2', 'kcal/mol/angstrom^2'),
            ],
            'topology':
                {
                    'type': 'bond',
                    'n_atoms': 2,
                    'symmetry': 'like_bond',
                    'fill': 0,
                    'flip': 0
                }
        },
    'quartic_bond':
        {
            'equation': ['K2*(R-R0)^2 + K3*(R-R0)^3 + K4*(R-R0)^4'],
            'constants':
                [
                    ('R0', 'angstrom'),
                    ('K2', 'kcal/mol/angstrom^2'),
                    ('K3', 'kcal/mol/angstrom^3'),
                    ('K4', 'kcal/mol/angstrom^4'),
                ],
            'topology':
                {
                    'type': 'bond',
                    'n_atoms': 2,
                    'symmetry': 'like_bond',
                    'fill': 0,
                    'flip': 0
                }
        },
    'quadratic_angle':
        {
            'equation': ['K2*(Theta-Theta0)^2'],
            'constants': [
                ('Theta0', 'degree'),
                ('K2', 'kcal/mol/radian^2'),
            ],
            'topology':
                {
                    'type': 'angle',
                    'n_atoms': 3,
                    'symmetry': 'like_angle',
                    'fill': 0,
                    'flip': 0
                }
        },
    'quartic_angle':
        {
            'equation':
                [
                    (
                        'K2*(Theta-Theta0)^2 + K3*(Theta-Theta0)^3'
                        '+ K4*(Theta-Theta0)^4'
                    )
                ],
            'constants':
                [
                    ('Theta0', 'degree'),
                    ('K2', 'kcal/mol/radian^2'),
                    ('K3', 'kcal/mol/radian^3'),
                    ('K4', 'kcal/mol/radian^4'),
                ],
            'topology':
                {
                    'type': 'angle',
                    'n_atoms': 3,
                    'symmetry': 'like_angle',
                    'fill': 0,
                    'flip': 0
                }
        },
    'torsion_1':
        {
            'equation': ['KPhi * [1 + cos(n*Phi - Phi0)]'],
            'constants': [
                ('KPhi', 'kcal/mol'),
                ('n', ''),
                ('Phi0', 'degree'),
            ],
            'topology':
                {
                    'type': 'torsion',
                    'n_atoms': 4,
                    'symmetry': 'like_torsion',
                    'fill': 0,
                    'flip': 0
                }
        },
    'torsion_2':
        {
            'equation': ['KPhi * [1 + cos(n*Phi - Phi0)]'],
            'constants': [
                ('KPhi', 'kcal/mol'),
                ('n', ''),
                ('Phi0', 'degree'),
            ],
            'topology':
                {
                    'type': 'torsion',
                    'n_atoms': 4,
                    'symmetry': 'like_torsion',
                    'fill': 0,
                    'flip': 0
                }
        },
    'torsion_3':
        {
            'equation':
                [
                    (
                        'V1 * [1 + cos(Phi - Phi0_1)]'
                        ' + V2 * [1 + cos(2*Phi - Phi0_2)]'
                        ' + V3 * [1 + cos(3*Phi - Phi0_3)]'
                    )
                ],
            'constants':
                [
                    ('V1', 'kcal/mol'),
                    ('Phi0_1', 'degree'),
                    ('V2', 'kcal/mol'),
                    ('Phi0_2', 'degree'),
                    ('V3', 'kcal/mol'),
                    ('Phi0_3', 'degree'),
                ],
            'topology':
                {
                    'type': 'torsion',
                    'n_atoms': 4,
                    'symmetry': 'like_torsion',
                    'fill': 0,
                    'flip': 0
                }
        },
    'wilson_out_of_plane':
        {
            'equation': ['K*(Chi - Chi0)^2'],
            'constants': [
                ('K', 'kcal/mol/radian^2'),
                ('Chi0', 'degree'),
            ],
            'topology':
                {
                    'type': 'out-of-plane',
                    'n_atoms': 4,
                    'symmetry': 'like_oop',
                    'fill': 0,
                    'flip': 0
                }
        },
    'nonbond(9-6)':
        {
            'equation':
                [
                    'eps(ij) [2(r(ij)*/r(ij))**9 - 3(r(ij)*/r(ij))**6]',
                    'r(ij) = [(r(i)**6 + r(j)**6))/2]**(1/6)',
                    (
                        'eps(ij) = 2 * sqrt(eps(i) * eps(j)) * '
                        'r(i)^3 * r(j)^3/[r(i)^6 + r(j)^6]'
                    )
                ],
            'constants': [('rmin', 'angstrom'), ('eps', 'kcal/mol')],
            'topology':
                {
                    'form': 'rmin-eps',
                    'type': 'pair',
                    'subtype': 'LJ 6-9',
                    'n_atoms': 1,
                    'symmetry': 'none',
                    'fill': 0,
                    'flip': 0
                }
        },
    'nonbond(12-6)':
        {
            'equation':
                [
                    'E = 4 * eps * [(sigma/r)**12 - (sigma/r)**6]',
                    'E = eps * [(rmin/r)**12 - (rmin/r)**6]',
                    'E = A/r**12 - B/r**6',
                    'rmin = 2**1/6 * sigma ',
                    'sigma = rmin / 2**1/6',
                    'A = 4 * eps * sigma**12',
                    'B = 4 * eps * sigma**6',
                    'sigma = (A/B)**1/6',
                    'eps = B**2/(4*A)'
                ],
            'constants': [('sigma', 'angstrom'), ('eps', 'kcal/mol')],
            'topology':
                {
                    'form': 'sigma-eps',
                    'type': 'pair',
                    'subtype': 'LJ 12-6',
                    'n_atoms': 1,
                    'symmetry': 'none',
                    'fill': 0,
                    'flip': 0
                }
        },
    'bond-bond':
        {
            'equation': ["K*(R-R0)*(R'-R0')"],
            'constants': [('K', 'kcal/mol/angstrom^2')],
            'topology':
                {
                    'type': 'bond-bond',
                    'n_atoms': 3,
                    'symmetry': 'like_angle',
                    'fill': 0,
                    'flip': 0
                }
        },
    'bond-bond_1_3':
        {
            'equation': ["K*(R-R0)*(R'-R0')"],
            'constants': [('K', 'kcal/mol/angstrom^2')],
            'topology':
                {
                    'type': '1,3 bond-bond',
                    'n_atoms': 4,
                    'symmetry': 'like_torsion',
                    'fill': 0,
                    'flip': 0
                }
        },
    'bond-angle':
        {
            'equation': ["K*(R-R0)*(Theta-Theta0)"],
            'constants':
                [
                    ('K12', 'kcal/mol/angstrom/radian'),
                    ('K23', 'kcal/mol/angstrom/radian'),
                ],
            'topology':
                {
                    'type': 'bond-angle',
                    'n_atoms': 3,
                    'symmetry': 'like_angle',
                    'fill': 1,
                    'flip': 1
                }
        },
    'angle-angle':
        {
            'equation': ["K*(Theta-Theta0)*(Theta'-Theta0')"],
            'constants': [('K', 'kcal/mol/angstrom/radian')],
            'topology':
                {
                    'type': 'angle-angle',
                    'n_atoms': 4,
                    'symmetry': 'like_angle-angle',
                    'fill': 0,
                    'flip': 0
                }
        },
    'end_bond-torsion_3':
        {
            'equation':
                [
                    (
                        '(R_L - R0_L) * (V1_L * [1 + cos(Phi - Phi0_1)]'
                        ' + V2_L * [1 + cos(2*Phi - Phi0_2)]'
                        ' + V3_L * [1 + cos(3*Phi - Phi0_3)])'
                    ),
                    (
                        '(R_R - R0_R) * (V1_R * [1 + cos(Phi - Phi0_1)]'
                        ' + V2_R * [1 + cos(2*Phi - Phi0_2)]'
                        ' + V3_R * [1 + cos(3*Phi - Phi0_3)])'
                    )
                ],
            'constants':
                [
                    ('V1_L', 'kcal/mol'),
                    ('V2_L', 'kcal/mol'),
                    ('V3_L', 'kcal/mol'),
                    ('V1_R', 'kcal/mol'),
                    ('V2_R', 'kcal/mol'),
                    ('V3_R', 'kcal/mol'),
                ],
            'topology':
                {
                    'type': 'torsion-end bond',
                    'n_atoms': 4,
                    'symmetry': 'like_torsion',
                    'fill': 3,
                    'flip': 3
                }
        },
    'middle_bond-torsion_3':
        {
            'equation':
                [
                    (
                        '(R_M - R0_M) * (V1 * [1 + cos(Phi - Phi0_1)]'
                        ' + V2 * [1 + cos(2*Phi - Phi0_2)]'
                        ' + V3 * [1 + cos(3*Phi - Phi0_3)])'
                    )
                ],
            'constants':
                [
                    ('V1', 'kcal/mol'),
                    ('V2', 'kcal/mol'),
                    ('V3', 'kcal/mol'),
                ],
            'topology':
                {
                    'type': 'torsion-middle bond',
                    'n_atoms': 4,
                    'symmetry': 'like_torsion',
                    'fill': 0,
                    'flip': 0
                }
        },
    'angle-torsion_3':
        {
            'equation':
                [
                    (
                        '(Theta_L - Theta0_L)'
                        '* (V1_L * [1 + cos(Phi - Phi0_1)]'
                        ' + V2_L * [1 + cos(2*Phi - Phi0_2)]'
                        ' + V3_L * [1 + cos(3*Phi - Phi0_3)])'
                    ),
                    (
                        '(Theta_R - Theta0_R)'
                        ' * (V1_R * [1 + cos(Phi - Phi0_1)]'
                        ' + V2_R * [1 + cos(2*Phi - Phi0_2)]'
                        ' + V3_R * [1 + cos(3*Phi - Phi0_3)])'
                    )
                ],
            'constants':
                [
                    ('V1_L', 'kcal/mol'),
                    ('V2_L', 'kcal/mol'),
                    ('V3_L', 'kcal/mol'),
                    ('V1_R', 'kcal/mol'),
                    ('V2_R', 'kcal/mol'),
                    ('V3_R', 'kcal/mol'),
                ],
            'topology':
                {
                    'type': 'torsion-angle',
                    'n_atoms': 4,
                    'symmetry': 'like_torsion',
                    'fill': 3,
                    'flip': 3
                }
        },
    'angle-angle-torsion_1':
        {
            'equation':
                [
                    (
                        'K * (Theta_L - Theta0_L) * (Theta_R - Theta0_R) * '
                        '(Phi - Phi0_1)'
                    )
                ],
            'constants': [('K', 'kcal/mol/degree^2/degree')],
            'topology':
                {
                    'type': 'angle-torsion-angle',
                    'n_atoms': 4,
                    'symmetry': 'like_torsion',
                    'fill': 0,
                    'flip': 0
                }
        },
    'torsion-torsion_1':
        {
            'equation': ['K * cos(Phi_L) * cos(Phi_R)'],
            'constants': [('K', 'kcal/mol')],
            'topology':
                {
                    'type': 'torsion-torsion',
                    'n_atoms': 5,
                    'symmetry': 'like_torsion-torsion',
                    'fill': 0,
                    'flip': 0
                }
        },
}
# yapf: enable
