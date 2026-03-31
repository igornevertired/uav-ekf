from typing import List


class Solution:
    def productExceptSelf(self, nums: List[int]) -> List[int]:
        answer = []
        mapl = 1
        mapr = 1
        for i in range(len(nums)):
            answer.append(mapl)
            mapl *= nums[i]
        print(answer)
        for i in range(len(answer) - 1, -1, -1):
            answer[i] *= mapr
            mapr *= nums[i]
        print(answer)


sol = Solution()
sol.productExceptSelf(nums=[1,2,3,4])