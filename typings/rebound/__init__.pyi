from typing import Any, Union, Sequence, Optional, Iterator, overload

class Orbit:
    a: float
    e: float
    inc: float
    Omega: float
    omega: float
    M: float
    f: float
    P: float
    d: float
    v: float

class Particle:
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    ax: float
    ay: float
    az: float
    m: float
    r: float
    name: str
    index: int
    def orbit(self, primary: Optional["Particle"] = ...) -> Orbit: ...

class Particles:
    @overload
    def __getitem__(self, key: Union[int, str]) -> Particle: ...
    @overload
    def __getitem__(self, key: slice) -> Sequence[Particle]: ...
    def __getitem__(self, key: Union[int, str, slice]) -> Any: ...
    def __setitem__(self, key: Union[int, str], value: Particle) -> None: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[Particle]: ...

class Simulation:
    integrator: str
    units: Sequence[str]
    G: float
    t: float
    dt: float
    N: int
    N_active: int
    N_var: int
    particles: Particles

    def __init__(self) -> None: ...
    def add(
        self,
        m: float = 0.0,
        r: float = 0.0,
        x: float = 0.0,
        y: float = 0.0,
        z: float = 0.0,
        vx: float = 0.0,
        vy: float = 0.0,
        vz: float = 0.0,
        a: float = 0.0,
        e: float = 0.0,
        inc: float = 0.0,
        Omega: float = 0.0,
        omega: float = 0.0,
        M: float = 0.0,
        primary: Optional[Particle] = None,
        name: str = "",
        **kwargs: Any
    ) -> None: ...
    def move_to_com(self) -> None: ...
    def move_to_hel(self) -> None: ...
    def integrate(self, time: float, exact_finish_time: int = 1) -> None: ...
    def energy(self) -> float: ...
