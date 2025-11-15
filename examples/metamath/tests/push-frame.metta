; MeTTa version of push-frame
; Original from mmverify-utils.metta

(= (adjust-state-additive $state $delta)
   (chain (get-state $state) $current
     (change-state! $state (+ $current $delta))))

(= (push-frame $fd)
  (adjust-state-additive $fd 1))

; Test 1: Push from 0 to 1
!(bind! &fd (new-state 0))
!(push-frame &fd)
!(assertEqual (get-state &fd) 1)

; Test 2: Push from 5 to 6
!(bind! &fd2 (new-state 5))
!(push-frame &fd2)
!(assertEqual (get-state &fd2) 6)
