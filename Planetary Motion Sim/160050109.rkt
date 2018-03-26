 #lang racket
(require "declarations.rkt")
(provide buildTree calcForces moveparticles)

(define (check-pos particle1 box)
   (let*[(cap (particle-posn particle1))]
        (cond[(and(>= (vec-x cap) (bbox-llx box)) (>= (vec-y cap) (bbox-lly box))
              (< (vec-x cap) (bbox-rux box))(< (vec-y cap) (bbox-ruy box))) #t]
             [else #f])))

(define (buildTree area config)
  (let*[(zeke (particle 0 (vec 0 0) (vec 0 0)))]
  (cond[(null? config) zeke]
       [(singleton config) (car config)]
       [else (gnode (netMass config area)
                    (COMass config)
                    (list (buildTree
                           (newArea1 area)(newConfig config '() (newArea1 area) 0))
                          (buildTree
                           (newArea2 area)(newConfig config '() (newArea2 area) 0))
                          (buildTree
                           (newArea3 area)(newConfig config '() (newArea3 area) 0))
                          (buildTree
                           (newArea4 area)(newConfig config '() (newArea4 area) 0))))])))
                          
(define (newArea1 area)
   (bbox (bbox-llx area) (* .5 (+ (bbox-lly area) (bbox-ruy area)))
              (* .5 (+ (bbox-llx area) (bbox-rux area))) (bbox-ruy area)))

(define (newArea2 area)
        (bbox (* .5 (+ (bbox-llx area) (bbox-rux area))) (* .5 (+ (bbox-lly area) (bbox-ruy area)))
              (bbox-rux area) (bbox-ruy area)))

(define (newArea3 area)
     (bbox (bbox-llx area) (bbox-lly area)
     (* .5 (+ (bbox-llx area) (bbox-rux area))) (* .5 (+ (bbox-lly area) (bbox-ruy area)))))

(define (newArea4 area)
  (bbox (* .5 (+ (bbox-llx area) (bbox-rux area))) (bbox-lly area)
              (bbox-rux area) (* .5 (+ (bbox-lly area) (bbox-ruy area)))))

(define (newConfig config config1 area n)
  (cond[(>= n (length config)) config1]
       [(check-pos (list-ref config n) area)
        (newConfig config (cons (list-ref config n) config1) area (+ n 1))]
       [else (newConfig config config1 area (+ n 1))]
       ))

(define (combMass list2 sum)
  (cond[(null? list2) sum]
       [(particle? list2) (particle-mass list2)]
       [else (combMass (cdr list2) (+ sum (particle-mass (car list2))))]))

         
(define (netMass list1 area)
  (let*[(newlist (newConfig list1 '() area 0))]
    (combMass newlist 0)))
 
(define (COMass config)
  (define (COMasshelp config)
    (cond[(null? config)(vec 0 0)]
         [(particle? config)
           (vec
            (* (particle-mass config) (vec-x (particle-posn config)))
            (* (particle-mass config) (vec-y (particle-posn config))))]
         [else (vec
                (foldr + 0 (map (lambda(k)(vec-x (COMasshelp k))) config))
                (foldr + 0 (map (lambda(k)(vec-y (COMasshelp k))) config)))]))
  (let*[(Mass (netMass config (bounding-box2 config)))]
    (if(zero? Mass) 0
    (vec (/ (vec-x (COMasshelp config)) Mass) (/ (vec-y (COMasshelp config)) Mass)))))

(define (bounding-box2 config)
  (if(particle? config) (bbox (- (vec-x (particle-posn config)) 1)
                              (- (vec-x (particle-posn config)) 1)
                              (+ 1 (vec-y (particle-posn config)))
                              (+ 1 (vec-y (particle-posn config))))
     (bounding-box config)))

(define (distance part1 gnode1)
  (sqrt(+ (sqr (- (vec-x (particle-posn part1)) (vec-x (gnode-posn gnode1))))
          (sqr (- (vec-y (particle-posn part1)) (vec-y (gnode-posn gnode1)))))))

(define (force particleOn gnodeBy)
  (let*[(term
         (/ (* g (particle-mass particleOn)(gnode-mass gnodeBy))
            (expt (distance particleOn gnodeBy) 3)))]
   (vec (* term (- (vec-x (gnode-posn gnodeBy)) (vec-x (particle-posn particleOn))))
        (* term (- (vec-y (gnode-posn gnodeBy)) (vec-y (particle-posn particleOn)))))))

(define (vec-add v1 v2)
  (vec (+ (vec-x v1) (vec-x v2))(+ (vec-y v1) (vec-y v2))))

(define (vec-mult v1 k)
  (vec (* k (vec-x v1)) (* k (vec-y v1))))

(define (calcForces area tree list)
   (cond[(null? list) '()]
        [else (append
                  (cons (netForce (car list) (vec 0 0) tree (- (bbox-rux area) (bbox-llx area))) '())
                  (calcForces area tree (cdr list)))]))
      
(define (netForce particle1 current tree side)
   (cond[(equal? particle1 tree)
            current]
        [(particle? tree)
         (vec-add current (force particle1 (gnode (particle-mass tree) (particle-posn tree) '())))]
        [(singleton (gnode-subtrees tree))
             (vec-add current (force particle1 (gnode (particle-mass (car tree)) (particle-posn (car tree)) '())))]
        [(< theta (/ (distance particle1 tree) side))
         (vec-add current (force particle1 tree))]
        [else (vec-add
                   (vec-add (netForce particle1 current (list-ref (gnode-subtrees tree) 0) (/ side 2))
                            (netForce particle1 current (list-ref (gnode-subtrees tree) 1) (/ side 2)))
                  (vec-add (netForce particle1 current (list-ref (gnode-subtrees tree) 2) (/ side 2))
                           (netForce particle1 current (list-ref (gnode-subtrees tree) 3) (/ side 2))))]))

(define (moveparticles particles forces)
   (cond[(null? particles) '()]
        [else(append (cons
               (particle
                 (particle-mass (car particles))
                 (newPosn (particle-posn (car particles)) (particle-velocity (car particles)) (accel (car particles) (car forces)) timeslice)
                 (newVel (particle-velocity (car particles)) (accel (car particles) (car forces)) timeslice))
              '())
           (moveparticles (cdr particles) (cdr forces)))]))


(define (newVel oldVel accel t)
  (vec-add oldVel (vec-mult accel t)))

(define (newPosn oldPosn oldVel accel t)
  (vec-add oldPosn (vec-add (vec-mult oldVel t) (vec-mult accel (* 0.5 t t)))))

(define (accel particle1 force)
  (vec-mult force (/ 1 (particle-mass particle1))))

